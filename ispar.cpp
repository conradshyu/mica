/*
 * ISPAR.CPP
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
 * last updated on April 26, 2005
*/

// support class implementation
#include "ispar.h"
#include "pthread.h"

// force PHP to return immediately
#define CLOSE_PHP   { fclose(stdin); fclose(stdout); fclose(stderr); };

const unsigned int nMaxTHREAD = 4;
const unsigned int nMaxBUFFER = 8192;

// define the structure for digest data
typedef struct
{
    vector<int> forward;    // list of forward fragments
    vector<int> reverse;    // list of reverse fragments
    int fshort, rshort;     // shortest forward and reverse fragments
    string accession;       // accession number
    string locus;           // locus
    string organism;        // organism name
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
}   // end of Timer(); linux implementation

#endif

/*
 * sort by the shortest forward fragment in the ascending order
*/
bool SortShortestForwardA(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.fshort < _b.fshort));
}   // end of SortShortestForwardA()

/*
 * sort by the shortest reverse fragment in the ascending order
*/
bool SortShortestReverseA(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.rshort < _b.rshort));
}   // end of SortShortestReverseA()

/*
 * sort by the shortest forward fragment in the descending order
*/
bool SortShortestForwardD(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.fshort > _b.fshort));
}   // end of SortShortestForwardD()

/*
 * sort by the shortest reverse fragment in the descending order
*/
bool SortShortestReverseD(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.rshort > _b.rshort));
}   // end of SortShortestReverseD()

/*
 * sort by the first forward fragment in the ascending order
*/
bool SortForwardFragmentA(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.forward.front() < _b.forward.front()));
}   // end of SortForwardFragmentA()

/*
 * sort by the first reverse fragment in the ascending order
*/
bool SortReverseFragmentA(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.reverse.front() < _b.reverse.front()));
}   // end of SortReverseFragmentA()

/*
 * sort by the first forward fragment in the descending order
*/
bool SortForwardFragmentD(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.forward.front() > _b.forward.front()));
}   // end of SortForwardFragmentD()

/*
 * sort by the first reverse fragment in the descending order
*/
bool SortReverseFragmentD(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.reverse.front() > _b.reverse.front()));
}   // end of SortReverseFragmentD()

/*
 * sort by the organism name in the ascending order
*/
bool SortOrganismA(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.organism < _b.organism));
}   // end of SortOrganismA()

/*
 * sort by the organism name in the descending order
*/
bool SortOrganismD(
    const stRECORD& _a, const stRECORD& _b)
{
    return((_a.organism > _b.organism));
}   // end of SortOrganismD()

/*
 * write the output in the plain text format, write all fragments
 * format: forward, reverse, ..., accession, locus, name
*/
bool WriteTXTa(
    list<stRECORD>& _lst,   // fragment data storage
    CmdParam&       _cmd)  // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".txt";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)     // file cannot be opened
    {
        return(false);
    }

    int fc = _cmd.EndonucleaseCount();

    ofs << "Query returned " << _lst.size() << " record(s)." << endl;
    ofs << "Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << endl;
    ofs << "Restriction Enzyme(s):";

    // print out the list of restriction enzymes
    for (int e = 0; e < _cmd.EndonucleaseCount(); ++e)
    {
        ofs << " " << _cmd.GetEndonuclease(e);
    }

    ofs << endl << endl;
    ofs << "Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer." << endl << endl;

    // print out the headings
    for (int c = 0; c < fc; ++c)
    {
        ofs.setf(ios::right);
        ofs << "Forward Reverse ";
    }

    ofs << "Accession Locus      Organism" << endl;

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        ofs.setf(ios::right);

        // write the fragment lengths into the file
        for (int f = 0; f < fc; ++f)
        {
            ofs << setw(7) << (*i).forward[f] << " ";
            ofs << setw(7) << (*i).reverse[f] << " ";
        }

        ofs.setf(ios::left);
        ofs << setw(9) << (*i).accession << " ";
        ofs << setw(10) << (*i).locus << " " << (*i).organism << endl;
    }

    ofs.close(); return(true);
}   // end of WriteTXTa(); write all fragments

/*
 * write the output in the plain text format, only the shortest fragments
 * format: forward, reverse, accession, locus, name
*/
bool WriteTXTs(
    list<stRECORD>& _lst,   // fragment data storage
    CmdParam&       _cmd)  // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".txt";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)     // file cannot be opened
    {
        return(false);
    }

    ofs << "Query returned " << _lst.size() << " record(s)." << endl;
    ofs << "Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << endl;
    ofs << "Restriction Enzyme(s):";

    // print out the list of restriction enzymes
    for (int e = 0; e < _cmd.EndonucleaseCount(); ++e)
    {
        ofs << " " << _cmd.GetEndonuclease(e);
    }

    ofs << endl << endl;
    ofs << "Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer." << endl << endl;

    ofs.setf(ios::right);
    ofs << "Forward Reverse Accession Locus      Organism" << endl;

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        // write the fragment lengths into the file
        ofs << setw(7) << (*i).fshort << " ";
        ofs << setw(7) << (*i).rshort << " ";

        ofs.setf(ios::left);
        ofs << setw(9) << (*i).accession << " ";
        ofs << setw(10) << (*i).locus << " " << (*i).organism << endl;
    }

    ofs.close(); return(true);
}   // end of WriteTXTs(); write the shortest fragments

/*
 * write the output in CSV format; Excel and other database programs readable
 * write all fragments
 * format: forward, reverse, ..., accession, locus, name
*/
bool WriteCSVa(
    list<stRECORD>& _lst,   // fragment data storage
    CmdParam&       _cmd)  // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".csv";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)     // file cannot be opened
    {
        return(false);
    }

    int fc = _cmd.EndonucleaseCount();

    ofs << "\"Query returned " << _lst.size() << " record(s).\"" << endl;
    ofs << "\"Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << "\"" << endl;
    ofs << "\"Restriction Enzyme(s):";

    // print out the list of restriction enzymes
    for (int e = 0; e < _cmd.EndonucleaseCount(); ++e)
    {
        ofs << " " << _cmd.GetEndonuclease(e);
    }

    ofs << "\"" << endl << endl;
    ofs << "\"Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer.\"" << endl << endl;

    // print out the headings
    for (int c = 0; c < fc; ++c)
    {
        ofs << "\"Forward\",\"Reverse\",";
    }

    ofs << "\"Accession\",\"Locus\",\"Organism\"" << endl;

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        // write the fragment lengths into the file
        for (int f = 0; f < fc; ++f)
        {
            ofs << (*i).forward[f] << "," << (*i).reverse[f] << ",";
        }

        ofs << "\"" << (*i).accession << "\",";       // accession number
        ofs << "\"" << (*i).locus << "\",";           // locus name
        ofs << "\"" << (*i).organism << "\"" << endl; // organism name
    }

    ofs.close(); return(true);
}   // end of WriteCSVa(); write all fragments

/*
 * write the output in CSV format; Excel and other database programs readable
 * write only the shortest fragments
 * format: forward, reverse, accession, locus, name
*/
bool WriteCSVs(
    list<stRECORD>& _lst,   // fragment data storage
    CmdParam&       _cmd)  // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".csv";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)     // file cannot be opened
    {
        return(false);
    }

    ofs << "\"Query returned " << _lst.size() << " record(s).\"" << endl;
    ofs << "\"Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << "\"" << endl;
    ofs << "\"Restriction Enzyme(s):";

    // print out the list of restriction enzymes
    for (int e = 0; e < _cmd.EndonucleaseCount(); ++e)
    {
        ofs << " " << _cmd.GetEndonuclease(e);
    }

    ofs << "\"" << endl << endl;
    ofs << "\"Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer.\"" << endl << endl;
    ofs << "\"Forward\",\"Reverse\",\"Accession\",\"Locus\",\"Organism\"" << endl;

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        // write the fragment lengths into the file
        ofs << (*i).fshort << "," << (*i).rshort << ",";
        ofs << "\"" << (*i).accession << "\",";       // accession number
        ofs << "\"" << (*i).locus << "\",";           // locus name
        ofs << "\"" << (*i).organism << "\"" << endl; // organism name
    }

    ofs.close(); return(true);
}   // end of WriteCSVs(); write the shortest fragments

/*
 * write the output in PAT format
 * species enzyme enzyme enzyme
*/
bool WritePAT(
    list<stRECORD>& _lst,       // fragment data storage
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".pat";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)
    {
        return(false);
    }

    int fc = _cmd.EndonucleaseCount();

    ofs << "Species\t";

    // print out the headings
    for (int c = 0; c < fc; ++c)
    {
        ofs << _cmd.GetEndonuclease(c) << "\t";
    }

    ofs << endl;

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        ofs << (*i).organism;

        // write the fragment lengths into the file
        for (int f = 0; f < fc; ++f)
        {
            ofs << "\t" << (*i).forward[f];
        }

        ofs << endl;
    }

    ofs.close(); return(true);
}   // end of WritePAT(); support format for PAT (phylogentic assignment tools)

/*
 * write the output in CSV format and use the PHP script to generate the web
 * interface for display, write all fragments
 * format: forward, reverse, ..., accession, locus, name
*/
bool WriteDATa(
    list<stRECORD>& _lst,       // fragment data storage
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".dat";             // add and extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)
    {
        return(false);
    }

    int fc = _cmd.EndonucleaseCount();

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        // write the fragment lengths into the file
        for (int f = 0; f < fc; ++f)
        {
            ofs << (*i).forward[f] << "," << (*i).reverse[f] << ",";
        }

        ofs << "\"" << (*i).accession << "\",";       // accession number
        ofs << "\"" << (*i).locus << "\",";           // locus name
        ofs << "\"" << (*i).organism << "\"" << endl; // organism name
    }

    ofs.close(); return(true);
}   // end of WriteDATa(); write all fragments for the php script

/*
 * write the output in CSV format and use the PHP script to generate the web
 * interface for display, write only the shortest fragments
 * format: forward, reverse, accession, locus, name
*/
bool WriteDATs(
    list<stRECORD>& _lst,       // fragment data storage
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".dat";             // add and extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)
    {
        return(false);
    }

    // iterate through the entire list
    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        // write the fragment lengths into the file
        ofs << (*i).fshort << "," << (*i).rshort << ",";
        ofs << "\"" << (*i).accession << "\",";       // accession number
        ofs << "\"" << (*i).locus << "\",";           // locus name
        ofs << "\"" << (*i).organism << "\"" << endl; // organism name
    }

    ofs.close(); return(true);
}   // end of WriteDATs(); write the shortest fragments for the php script

/*
 * draw the output in PHP format for web display
*/
bool WritePHP(
    list<stRECORD>& _lst,       // fragment data storage
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".php";             // add and extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)
    {
        return(false);
    }

    int fc = (_cmd.OutputAll()) ? _cmd.EndonucleaseCount() : 1;
    fc = (fc * 2) + 3;

    ofs << "<?php" << endl;
    ofs << "  require \"shared.data.inc\";" << endl;
    ofs << "  DrawHeader(\"MiCA: Virtual Digest (ISPaR) Output\");" << endl;
    ofs << "  $page = $HTTP_GET_VARS[\'page\'];" << endl;
    ofs << "  DrawISPAR(" << fc << ", " << _cmd.GetFilename() << ", $page);" << endl;
    ofs << "  DrawFooter();" << endl << "?>" << endl;
    ofs.close(); return(true);
}   // end of WritePHP()

// globally accessible classes for multithreading
CmdParam cmd; SeqDB rdp;
list<stRECORD> data;
pthread_mutex_t mtxLockDBMS;    // critical region lock for database
pthread_mutex_t mtxLockITEM;    // critical region lock for records

/*
 * the prodcution trflp function
*/
void* DoDigest(void*)
{
    tRFLP rflp(cmd);
    stRECORD item; bool run = false;

    do
    {
        pthread_mutex_lock(&mtxLockDBMS);
        // ** enter the critical section for database
        run = rdp.NextRecord();     // retrive a sequence from the database

        if (run)
        {
            rflp.SetStrand(rdp.GetOrigin());
            item.locus = rdp.GetLocus(), item.organism = rdp.GetOrganism();
            item.accession = rdp.GetAccession();
        }   // set the sequence for the search
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLockDBMS);

        if (!run || !rflp.Delimit())
        {
            continue;
        }   // skip if there is no more sequences or amplification fails

        rflp.Digest();
        rflp.GetFragment(item.forward, item.reverse);     // all fragments
        rflp.GetFragment(item.fshort, item.rshort);       // shortest fragments

        pthread_mutex_lock(&mtxLockITEM);
        // ** enter the critical section for records
        data.push_back(item);
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLockITEM);
    } while (run);

    return(NULL);
}   // end of DoDigest(); production function for trflp

/*
 * to compile, type:
 * g++ cmdparam.cpp seqdb.cpp bitvector.cpp ispar.cpp -o ispar.exe
*/
int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "insufficient number of parameters" << endl;
        return(1);
    }

    cmd.OpenFile(argv[1]);              // open the parameter file
    rdp.OpenFile(cmd.GetDatabase());      // open the sequence database
    data.clear();

    pthread_mutex_init(&mtxLockDBMS, NULL);   // initialize the lock for database
    pthread_mutex_init(&mtxLockITEM, NULL);   // initialize the lock for record
    vector<pthread_t> pts(nMaxTHREAD, 0);     // a vector for pthread

    for (unsigned int i = 0; i < pts.size(); ++i)
    {
        pthread_create(&pts[i], NULL, &DoDigest, NULL);
    }   // first, spawn the threads

    for (unsigned int j = 0; j < pts.size(); ++j)
    {
        pthread_join(pts[j], NULL);
    }   // wait for all threads to complete

    bool (*SortOption[10])(const stRECORD&, const stRECORD&) =
    {
        SortForwardFragmentA,   // forward fragments in ascending order
        SortReverseFragmentA,   // reverse fragments in ascending order
        SortShortestForwardA,   // shortest forward fragment in ascending order
        SortShortestReverseA,   // shortest reverse fragment in ascending order
        SortOrganismA,          // organism name in ascending order
        SortForwardFragmentD,   // forward fragments in descending order
        SortReverseFragmentD,   // reverse fragments in descending order
        SortShortestForwardD,   // shortest forward fragment in descending order
        SortShortestReverseD,   // shortest reverse fragment in descending order
        SortOrganismD           // organism name in descending order
    };  // nasty function pointers

    bool (*WriteOutput[2][5])(list<stRECORD>&, CmdParam&) =
    {
        { WriteTXTa, WriteCSVa, WritePAT, WriteDATa, WritePHP },
        { WriteTXTs, WriteCSVs, WritePAT, WriteDATs, WritePHP }
    };

    data.sort(*SortOption[cmd.SortOption()]);     // sort the output data

    // write the query results in different formats
    WriteOutput[static_cast<int>(cmd.OutputShort())][0](data, cmd);
    WriteOutput[static_cast<int>(cmd.OutputShort())][1](data, cmd);
    WriteOutput[static_cast<int>(cmd.OutputShort())][2](data, cmd);
    WriteOutput[static_cast<int>(cmd.OutputShort())][3](data, cmd);
    WriteOutput[static_cast<int>(cmd.OutputShort())][4](data, cmd);

    pthread_exit(NULL);
    pthread_mutex_destroy(&mtxLockITEM);
    pthread_mutex_destroy(&mtxLockDBMS);
    return(0);
}   // end of main()
