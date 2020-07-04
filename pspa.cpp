/*
 * PSPA.CPP
 *
 * primer sequence prevalence analysis
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
 * last updated on April 28, 2005
 * last revised on July 10, 2007
*/

// support class implementation
#include "pspa.h"
#include "pthread.h"

// force PHP to return immediately
#define CLOSE_PHP   { fclose(stdin); fclose(stdout); fclose(stderr); };

const unsigned int nMaxTHREAD = 4;
const unsigned int nMaxBUFFER = 2048;

// define the structure for digest data
typedef struct
{
    string forward;         // forward primer sequence
    string reverse;         // reverse primer sequence
    int forward_match;      // number of matches for forward primer
    int reverse_match;      // number of matches for reverse primer
    int primer_match;       // number of matches for both primers
} stRECORD;

/*
 * non-class implementations
*/
/*
 * sort by the forward primer in the ascending order
*/
bool SortForwardPrimerA(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.forward < _b.forward);
}   // end of SortForwardPrimerA()

/*
 * sort by the reverse primer in the ascending order
*/
bool SortReversePrimerA(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.reverse < _b.reverse);
}   // end of SortReversePrimerA()

/*
 * sort by the number of both primers that has been found in the ascending order
*/
bool SortPrimerMatchA(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.primer_match < _b.primer_match);
}   // end of SortPrimerMatchA()

/*
 * sort by the forward primer in the descending order
*/
bool SortForwardPrimerD(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.forward > _b.forward);
}   // end of SortForwardPrimerD()

/*
 * sort by the reverse primer in the descending order
*/
bool SortReversePrimerD(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.reverse > _b.reverse);
}   // end of SortReversePrimerD()

/*
 * sort by the number of both primers that has been found in the descending order
*/
bool SortPrimerMatchD(
    const stRECORD& _a, const stRECORD& _b)
{
    return(_a.primer_match > _b.primer_match);
}   // end of SortPrimerMatchD()

/*
 * write the output in the plain text format
 * format: forward primer, reverse primer, matches
*/
bool WriteTXT(
    list<stRECORD>& _lst,   // fragment data storage
    CmdParam&       _cmd)  // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".txt";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    ofs << "Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer." << endl << endl;
    ofs << "Forward Primer            Forward Matches Reverse Primer            ";
    ofs <<" Reverse Matches Both Matches" << endl;

    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        sprintf(buffer, "%25s %15d %25s %15d %12d",
            (*i).forward.c_str(),     // forward primer sequence
            (*i).forward_match,       // number of forward primer amplified sequences
            (*i).reverse.c_str(),     // reverse primer sequence
            (*i).reverse_match,       // number of reverse primer amplified sequences
            (*i).primer_match);      // number of sequences amplified by both primers
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteTXT()

/*
 * write the output in CSV format; Excel and other database programs readable
 * write all fragments
*/
bool WriteCSV(
    list<stRECORD>& _lst,   // fragment data storage
    CmdParam&       _cmd)  // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".csv";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    ofs << "\"Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer.\"" << endl << endl;

    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        ofs << "\"" << (*i).forward << "\"," << (*i).forward_match;
        ofs << ",\"" << (*i).reverse << "\"," << (*i).reverse_match;
        ofs << "," << (*i).primer_match << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteCSV()

/*
 * write the output in CSV format and use the PHP script to generate the web
 * interface for display, write all fragments
*/
bool WriteDAT(
    list<stRECORD>& _lst,       // fragment data storage
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".dat";             // add and extension
    ofstream ofs(name.c_str(), ios::trunc);

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    for (list<stRECORD>::iterator i = _lst.begin(); !(i == _lst.end()); ++i)
    {
        ofs << "\"" << (*i).forward << "\"," << (*i).forward_match;
        ofs << ",\"" << (*i).reverse << "\"," << (*i).reverse_match;
        ofs << "," << (*i).primer_match << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteDAT()

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
    }   // file cannot be opened or created successfully

    ofs << "<?php" << endl;
    ofs << "  require \"shared.data.inc\";" << endl;
    ofs << "  DrawHeader(\"MiCA: Primer Sequence Prevalence Anlysis Output\");" << endl;
    ofs << "  DrawPSPA(" << _cmd.GetFilename() << ");" << endl;
    ofs << "  DrawFooter();" << endl << "?>" << endl;

    ofs.close(); return(true);
}   // end of WritePHP()

// globally accessible classes for multithreading
CmdParam cmd; SeqDB rdp;
list<stRECORD> data;
pthread_mutex_t mtxLockDBMS;    // critical region lock for database
pthread_mutex_t mtxLockITEM;    // critical region lock for records

/*
 * the prodcution pspa function
*/
void* DoPSPA(void*)
{
    cPSPA pspa(cmd);          // initialize the class
    vector<bool> forward, reverse;
    unsigned int idx;
    bool run = false;

    do
    {
        pthread_mutex_lock(&mtxLockDBMS);
        // ** enter the critical section for database
        run = rdp.NextRecord();     // retrive a sequence from the database

        if (run)
        {
            pspa.SetStrand(rdp.GetOrigin());
        }   // set the sequence for the search
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLockDBMS);

        if (!run)
        {
            continue;
        }   // skip if there is no more sequences

        forward.clear(); reverse.clear();
        pspa.Delimit(forward, reverse);   // perform the search on all primers
        idx = 0;

        pthread_mutex_lock(&mtxLockITEM);
        // ** enter the critical section for records
        for (list<stRECORD>::iterator k = data.begin(); !(k == data.end()); ++k, ++idx)
        {
            (*k).forward_match += static_cast<int>(forward[idx]);
            (*k).reverse_match += static_cast<int>(reverse[idx]);

            if (forward[idx] && reverse[idx])
            {
                (*k).primer_match += 1;
            }   // calculate the number of simultaneous matches
        }   // record the digestions
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLockITEM);
    } while (run);

    return(NULL);
}   // end of DoPSPA()

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
    data.clear(); stRECORD item;

    for (unsigned int f = 0; f < cmd.ForwardPrimerCount(); ++f)
    {
        for (unsigned int r = 0; r < cmd.ReversePrimerCount(); ++r)
        {
            item.forward = cmd.GetForwardPrimer(f);
            item.reverse = cmd.GetReversePrimer(r);
            item.forward_match = item.reverse_match = item.primer_match = 0;
            data.push_back(item);

#ifdef _VERBOSE
            cout << "forward primer: " << item.forward << endl;
            cout << "reverse primer: " << item.reverse << endl;
#endif  // _VERBOSE
        }
    }   // construct and initialize the records

    pthread_mutex_init(&mtxLockDBMS, NULL);   // initialize the lock for database
    pthread_mutex_init(&mtxLockITEM, NULL);   // initialize the lock for record
    vector<pthread_t> pts(nMaxTHREAD, 0);     // a vector for pthread

    for (unsigned int i = 0; i < pts.size(); ++i)
    {
        pthread_create(&pts[i], NULL, &DoPSPA, NULL);
    }   // first, spawn the threads

    for (unsigned int j = 0; j < pts.size(); ++j)
    {
        pthread_join(pts[j], NULL);
    }   // wait for all threads to complete

#ifdef _VERBOSE
    for (list<stRECORD>::iterator d = data.begin(); !(d == data.end()); ++d)
    {
        cout << (*d).forward << ", " << (*d).forward_match << ", ";
        cout << (*d).reverse << ", " << (*d).reverse_match << ", ";
        cout << (*d).primer_match << endl;
    }   // print out the content of list
#endif  // _VERBOSE

    bool (*SortOption[6])(const stRECORD&, const stRECORD&) =
    {
        SortForwardPrimerA, // sort by forward primers in ascending order
        SortReversePrimerA, // sort by reverse primers in ascending order
        SortPrimerMatchA,   // sort by number of both primers matches in ascending order
        SortForwardPrimerD, // sort by reverse primers in descending order
        SortReversePrimerD, // sort by reverse primers in descending order
        SortPrimerMatchD    // sort by number of primer matches in descending order
    };  // nasty function pointers

    data.sort(*SortOption[cmd.SortOption()]);   // sort the output data

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
