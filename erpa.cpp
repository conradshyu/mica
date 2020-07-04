/*
 * ERPA.CPP
 *
 * enzyme resolving power analysis
 *
 * Written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83844
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * this program performs terminal restriction fragment length polymorphism (T-RFLP)
 * analysis on 16S and 18S ribosomal RNA genes
 *
 * All rights reserved. Copyrights (R) 2005.
 * last updated on May 5, 2005
 * last revised on July 10, 2007
*/

// support class implementation
#include "erpa.h"
#include "pthread.h"

// force PHP to return immediately
//#define _VERBOSE
#define CLOSE_PHP       { fclose(stdin); fclose(stdout); fclose(stderr); };

const int nMaxTHREAD    = 4;
const int nMaxBUFFER    = 8192;
const int nMaxFRAGMENT  = 2000;

typedef struct
{
    vector<int> forward;    // histogram for forward fragments
    vector<int> reverse;    // histogram for reverse fragments
    string site;            // restriction site
} stHISTOGRAM;

// define the structure for digest data
typedef struct
{
    vector<int> forward;    // list of forward fragments
    vector<int> reverse;    // list of reverse fragments
    double forward_mean;    // average forward fragment sizes
    double forward_stdev;   // standard deviation for forward fragments
    double reverse_mean;    // average reverse fragment sizes
    double reverse_stdev;   // standard deviation for reverse fragments
    int reverse_unique;     // number of unique reverse fragments
    int forward_unique;     // number of unique forward fragments
    int success;            // number of successful cuts
    string site;            // sequence of the restriction enzyme
} stRECORD;

/*
 * non-class implementations
*/
#ifdef _WIN32

#include <windows.h>

// for Windows system; minimum operating system: windows 95, windows NT 3.1
double Timer()
{
    LARGE_INTEGER freq, timer;

    /*
     * the QueryPerformanceFrequency function retrieves the frequency of the high-
     * resolution performance counter, if one exists. The frequency cannot change
     * while the system is running.
    */
    QueryPerformanceFrequency(&freq);

    /*
     * the QueryPerformanceCounter function retrieves the current value of the
     * high-resolution performance counter.
    */
    QueryPerformanceCounter(&timer);

    return(static_cast<double>(timer.QuadPart / freq.QuadPart));
}   // end of Timer(); windows implementation

#else

// for Unix/Linux systems
#include <sys/time.h>

double Timer()
{
    struct timeval tm;

    /*
     * the gettimeofday function gets the current time, expressed as seconds and
     * microseconds since 00:00 Coordinated Universal Time (UTC), January 1, 1970.
    */
    gettimeofday(&tm, 0);

    return(static_cast<double>(tm.tv_sec + tm.tv_usec / 1000000.0));
}   // end of Timer(); linux or unix implementation

#endif

/*
 * sort by restriction enzymes in ascending order
*/
bool SortEndonucleaseA(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.site < _b.site);
}   // end of SortEndonucleaseA()

/*
 * sort by the number of unique forward fragments in ascending order
*/
bool SortForwardUniqueA(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.forward_unique < _b.forward_unique);
}   // end of SortForwardUniqueA()

/*
 * sort by the number of unique reverse fragments in descending order
*/
bool SortReverseUniqueA(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.reverse_unique < _b.reverse_unique);
}   // end of SortReverseUniqueA()

/*
 * sort by the restriction enzymes in the descending order
*/
bool SortEndonucleaseD(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.site > _b.site);
}   // end of SortEndonucleaseD()

/*
 * sort by the number of unique forward fragments in the descending order
*/
bool SortForwardUniqueD(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.forward_unique > _b.forward_unique);
}   // end of SortForwardUniqueD()

/*
 * sort by the number of unique reverse fragments in the descending order
*/
bool SortReverseUniqueD(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.reverse_unique > _b.reverse_unique);
}   // end of SortReverseUniqueD()

/*
 * write the output in the plain text format
 *
 * format:
 * Restrict Site Total Hits Full Length 5'unique 5'average 5'stddev 3'unique 3'average 3'stddev
 *    TC^NNGA    1234567890 12345678901 12345678 123456789 12345678 12345678 123456789 12345678
*/
bool WriteTXT(
    list<stRECORD>& _lst,
    CmdParam&       _cmd)
{
    string name = _cmd.GetFilename();
    name += ".txt";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    ofs << "Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << endl << endl;
    ofs << "Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer." << endl << endl;
    ofs << "Restrict Site Total Hits Full Length ";
    ofs << "5\'Unique 5\'Average 5\'StdDev 3\'Unique 3\'Average 3\'StdDev" << endl;

    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        sprintf(buffer, "%13s %10d %11d %8d %9.3f %8.3f %8d %9.3f %8.3f",
            (*i).site.c_str(), static_cast<int>((*i).forward.size()), (*i).success,
            (*i).forward_unique, (*i).forward_mean, (*i).forward_stdev,
            (*i).reverse_unique, (*i).reverse_mean, (*i).reverse_stdev);
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteTXT()

/*
 * write the output in the CSV (comma separated value) format
 *
 * format:
 * Restrict Site Total Hits Full Length 5'unique 5'average 5'stddev 3'unique 3'average 3'stddev
 *    TC^NNGA    1234567890 12345678901 12345678 123456789 12345678 12345678 123456789 12345678
*/
bool WriteCSV(
    list<stRECORD>& _lst,
    CmdParam&       _cmd)
{
    string name = _cmd.GetFilename();
    name += ".csv";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];  // string buffer for formating

    if (!ofs)     // file cannot be opened
    {
        return(false);
    }

    ofs << "\"Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << "\"" << endl << endl;
    ofs << "\"Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer.\"" << endl << endl;
    ofs << "\"Restrict Site\",\"Total Hits\",\"Full Length\",\"5\'Unique\",";
    ofs << "\"5\'Average\",\"5\'StdDev\",\"3\'Unique\",\"3\'Average\",\"3\'StdDev\"" << endl;

    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        sprintf(buffer, "\"%s\",%d,%d,%d,%.3f,%.3f,%d,%.3f,%.3f",
            (*i).site.c_str(), static_cast<int>((*i).forward.size()), (*i).success,
            (*i).forward_unique, (*i).forward_mean, (*i).forward_stdev,
            (*i).reverse_unique, (*i).reverse_mean, (*i).reverse_stdev);
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteCSV()

/*
 * write the output in the CSV (comma separated value) format for the php script
*/
bool WriteDAT(
    list<stRECORD>& _lst,
    CmdParam&       _cmd)
{
    string name = _cmd.GetFilename();
    name += ".dat";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        sprintf(buffer, "\"%s\",%d,%d,%d,%.3f,%.3f,%d,%.3f,%.3f",
            (*i).site.c_str(), static_cast<int>((*i).forward.size()), (*i).success,
            (*i).forward_unique, (*i).forward_mean, (*i).forward_stdev,
            (*i).reverse_unique, (*i).reverse_mean, (*i).reverse_stdev);
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteDAT()

/*
 * write the php script for web display
*/
bool WritePHP(
    list<stRECORD>& _lst,
    CmdParam&       _cmd)
{
    string name = _cmd.GetFilename();
    name += ".php";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)     // file cannot be opened
    {
        return(false);
    }

    ofs << "<?php" << endl;
    ofs << "  require \"shared.data.inc\";" << endl;
    ofs << "  DrawHeader(\"MiCA: Enzyme Resolving Power Analysis Output\");" << endl;
    ofs << "  DrawERPA(" << _cmd.GetFilename() << ");" << endl;
    ofs << "  DrawFooter();" << endl;
    ofs << "?>" << endl;

    ofs.close(); return(true);
}   // end of WritePHP()

/*
 * calculate the number of unique elements in the vector
*/
int Unique(
    vector<int>& _vtr)
{
    vector<bool> element(nMaxFRAGMENT, false);
    int unique = 0;

    for (unsigned int i = 0; i < _vtr.size(); ++i)
    {
        if (!(element[_vtr[i]]))
        {
            element[_vtr[i]] = true; ++unique;
        }
    }   // iterate through the entire vector

    return(unique);
}   // end of Unique()

/*
 * calculate the histogram of the fragments
*/
bool Histogram(
    list<stRECORD>& _rec)
{
    for (list<stRECORD>::iterator s = _rec.begin(); !(s == _rec.end()); ++s)
    {
        cout << "\"" << (*s).site << "\",,";
    }   // iterate through the entire list of records and calculate the histogram

    cout << endl;

    for (unsigned int v = 0; v < (_rec.front()).forward.size() ; ++v)
    {
        for (list<stRECORD>::iterator r = _rec.begin(); !(r == _rec.end()); ++r)
        {
            cout << (*r).forward[v] << "," << (*r).reverse[v] << ",";
        }

        cout << endl;
    }

    return(true);
}   // end of Histogram()

/*
 * calculate the mean and standard deviation; this function will calculate the
 * unbiased estimate of the sample variance
*/
bool Statistics(
    vector<int>& _vtr,      // vector of numbers
    double& _mean,          // mean of the numbers
    double& _stdev)        // unbiased estimate of the sample variance
{
    if (!(_vtr.size()))
    {
        return(false);
    }

    _mean = _stdev = 0.0;

    for (unsigned int i = 0; i < _vtr.size(); ++i)
    {
        _mean += _vtr[i];
    }   // first calculate the mean

    _mean /= _vtr.size();

    for (unsigned int j = 0; j < _vtr.size(); ++j)
    {
        _stdev += pow((_mean - _vtr[j]), 2.0);
    }   // accumulate the differences

    _stdev = sqrt(_stdev / (_vtr.size() - 1));  return(true);
}   // end of Statistics()

// globally accessible classes for multithreading
CmdParam cmd; SeqDB rdp;
list<stRECORD> data;
pthread_mutex_t mtxLockDBMS;    // critical region lock for database
pthread_mutex_t mtxLockITEM;    // critical region lock for records

/*
 * perform the enzyme resolving power analysis
*/
void* DoERPA(void*)
{
    cERPA rflp(cmd);      // instantiate the class
    int forward, reverse; bool run = false;
    int success;

    do
    {
        pthread_mutex_lock(&mtxLockDBMS);
        // ** enter the critical section for database
        run = rdp.NextRecord();     // retrive a sequence from the database

        if (run)
        {
            rflp.SetStrand(rdp.GetOrigin());
        }   // set the sequence for the search
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLockDBMS);

        if (!run || !rflp.Delimit())
        {
            continue;
        }   // skip if there is no more sequences or amplification failed

        for (list<stRECORD>::iterator k = data.begin(); !(k == data.end()); ++k)
        {
            // accumulate the number of successful cuts; true = 1; false = 0
            success = static_cast<int>(rflp.Digest((*k).site));

            pthread_mutex_lock(&mtxLockITEM);
            // ** enter the critical section for records
            (*k).success += success;
            rflp.GetFragment(forward, reverse);       // get the fragment sizes
            (*k).forward.push_back(forward);        // store the forward fragment
            (*k).reverse.push_back(reverse);        // store the reverse fragment
            // ** leave the critical section for database
            pthread_mutex_unlock(&mtxLockITEM);
        }   // iterate through the entire list of restriction enzymes
    } while (run);

    return(NULL);
}   // end of DoERPA()

/*
 * to compile, type:
 * g++ cmdparam.cpp seqdb.cpp bitvector.cpp ispar.cpp -o ispar.exe
*/
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "insufficient number of parameters" << std::endl;
        return(1);
    }

    cmd.OpenFile(argv[1]);              // open the parameter file
    rdp.OpenFile(cmd.GetDatabase());      // open the sequence database
    data.clear(); stRECORD item;

    for (unsigned int e = 0; e < cmd.EndonucleaseCount(); ++e)
    {
        item.forward.clear(); item.reverse.clear();
        item.forward_mean = item.forward_stdev = 0.0;
        item.reverse_mean = item.reverse_stdev = 0.0;
        item.forward_unique = item.reverse_unique = 0;
        item.success = 0; item.site = cmd.GetEndonuclease(e);
        data.push_back(item);
    }   // initialize and record the restriction site

    pthread_mutex_init(&mtxLockDBMS, NULL);   // initialize the lock for database
    pthread_mutex_init(&mtxLockITEM, NULL);   // initialize the lock for record
    vector<pthread_t> pts(nMaxTHREAD, 0);     // a vector for pthread

    for (unsigned int i = 0; i < pts.size(); ++i)
    {
        pthread_create(&pts[i], NULL, &DoERPA, NULL);
    }   // first, spawn the threads

    for (unsigned int j = 0; j < pts.size(); ++j)
    {
        pthread_join(pts[j], NULL);
    }   // wait for all threads to complete

#ifdef _FRAGMENTS   // output all fragements for the creations of histograms
    Histogram(data);
#endif  // _FRAGMENTS

    for (list<stRECORD>::iterator s = data.begin(); !(s == data.end()); ++s)
    {
        Statistics((*s).forward, (*s).forward_mean, (*s).forward_stdev);
        Statistics((*s).reverse, (*s).reverse_mean, (*s).reverse_stdev);
        (*s).forward_unique = Unique((*s).forward);
        (*s).reverse_unique = Unique((*s).reverse);

#ifdef _VERBOSE
        cout << (*s).site << ", " << (*s).success << ", " << (*s).forward.size();
        cout << ", " << (*s).forward_mean << ", " << (*s).forward_stdev;
        cout << ", " << (*s).forward_unique << ", " << (*s).reverse_mean;
        cout << ", " << (*s).reverse_stdev << ", " << (*s).reverse_unique << endl;
#endif  // _VERBOSE
    }   // iterate through the entire list of records and calculate the average lengths


    bool (*SortOption[6])(const stRECORD&, const stRECORD&) =
    {
        SortEndonucleaseA,      // sort by restriction enzymes in ascending order
        SortForwardUniqueA,     // sort by forward fragments in ascending order
        SortReverseUniqueA,     // sort by reverse fragments in descending order
        SortEndonucleaseD,      // sort by restriction enzymes in descending order
        SortForwardUniqueD,     // sort by forward fragments in ascending order
        SortReverseUniqueD      // sort by reverse fragments in descending order
    };  // nasty function pointers

    data.sort(*SortOption[cmd.SortOption()]);     // sort the output data

    /*
     * write the output in various formats; explicitly signal the compiler that these
     * funcation calls do not require any specific order. the compiler is free to
     * rearrange for processor dispatch and parallel processing.
    */
    WriteTXT(data, cmd), WriteCSV(data, cmd), WriteDAT(data, cmd), WritePHP(data, cmd);

    pthread_exit(NULL);
    pthread_mutex_destroy(&mtxLockITEM);
    pthread_mutex_destroy(&mtxLockDBMS);
    return(0);
}   // end of main()
