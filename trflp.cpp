/*
 * TRFLP.CPP
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
 * revised on November 15, 2005
 * major revision on the calculation of abundance on July 9, 2007
 * last revised on July 10, 2007
 *
 * Note: added support for using only forward or reverse fragments for the
 * identification of species composition
*/

// support class implementation
#include "trflp.h"
#include "pthread.h"

// force PHP to return immediately
#define CLOSE_PHP       { fclose(stdin); fclose(stdout); fclose(stderr); };

const unsigned int nMaxTHREAD = 4;

/*
 * sort by the species abundance in ascending order
*/
bool SortBiomassA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.biomass > _b.biomass));
}   // end of SortBiomassA()

/*
 * sort by the species abundance in descending order
*/
bool SortBiomassD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.biomass < _b.biomass));
}   // end of SortBiomassD

/*
 * sort by the species name in ascending order
*/
bool SortOrganismA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.organism > _b.organism));
}   // end of SortOrganismA()

/*
 * sort by the species name is descending order
*/
bool SortOrganismD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.organism < _b.organism));
}   // end of SortOrganismD()

/*
 * sort by the sample forward fragment size in ascending order
*/
bool SortForwardA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.fobserve > _b.fobserve));
}   // end of SortForwardA()

/*
 * sort by the sample forward fragment size in descending order
*/
bool SortForwardD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.fobserve < _b.fobserve));
}   // end of SortForwardD()

/*
 * sort by the sample reverse fragment size in ascending order
*/
bool SortReverseA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.robserve > _b.robserve));
}   // end of SortReverseA()

/*
 * sort by the sample reverse fragment size in descending order
*/
bool SortReverseD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.robserve < _b.robserve));
}   // end of SortReverseD()

/*
 * write the output in the plain text format
 * format:
 * predicted forward and revesre fragments
 * observed (sample) forward and reverse fragments
 * normalized relative abundance
 * organism name
 *
 * last updated on July 14, 2006
*/
bool WriteDAT(
    list<stNICHE>&  _niche,     // plausible community profile based on trflp data
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".dat";         // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    for (list<stNICHE>::iterator i = _niche.begin(); !(i == _niche.end()); ++i)
    {
        sprintf(buffer, "%.0f,%.0f,%.2f,%.2f,%.6f,\"%s\"",
            (*i).fpredict,            // predicted forward fragment
            (*i).rpredict,            // predicted reverse fragment
            (*i).fobserve,            // observed (sample) forward fragment
            (*i).robserve,            // observed (sample) reverse fragment
            (*i).biomass,             // normalized relative abundance
            (*i).organism.c_str());  // name of the species
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteDAT()

/*
 * write the output in the plain text format, only the shortest fragments
 * format:
 * predicted forward and revesre fragments
 * observed (sample) forward and reverse fragments
 * normalized relative abundance
 * organism name
 *
 * last updated on July 14, 2006
 * last revised on July 10, 2007
*/
bool WriteTXT(
    list<stNICHE>&  _niche,     // fragment data storage
    CmdParam&       _cmd)      // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".txt";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    ofs << "Query returned " << _niche.size() << " record(s)." << endl;
    ofs << "Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << endl;
    ofs << "Restriction Enzyme: " << _cmd.GetEndonuclease(0) << endl;
    ofs << endl << "Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer." << endl << endl;

    ofs.setf(ios::right);
    ofs << "Forward Reverse Forward Reverse Abundance Name " << endl;

    for (list<stNICHE>::iterator i = _niche.begin(); !(i == _niche.end()); ++i)
    {
        sprintf(buffer, "%7.0f %7.0f %7.2f %7.2f %9.6f %s",
            (*i).fpredict,            // predicted forward fragment
            (*i).rpredict,            // predicted reverse fragment
            (*i).fobserve,            // observed (sample) forward fragment
            (*i).robserve,            // observed (sample) reverse fragment
            (*i).biomass,             // normalized relative abundance
            (*i).organism.c_str());  // name of the species
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteTXT()

/*
 * write the output in CSV format; Excel and other database programs readable
 * format:
 * predicted forward and revesre fragments
 * observed (sample) forward and reverse fragments
 * normalized relative abundance
 * organism name
 *
 * last updated on July 14, 2006
*/
bool WriteCSV(
    list<stNICHE>& _niche,      // fragment data storage
    CmdParam&      _cmd)       // command-line parameters
{
    string name = _cmd.GetFilename();
    name += ".csv";             // add an extension
    ofstream ofs(name.c_str(), ios::trunc);
    char buffer[nMaxBUFFER];

    if (!ofs)
    {
        return(false);
    }   // file cannot be opened or created successfully

    ofs << "\"Query returned " << _niche.size() << " record(s).\"" << endl;
    ofs << "\"Forward Primer: " << _cmd.GetForwardPrimer(0) << ", ";
    ofs << "Reverse Primer: " << _cmd.GetReversePrimer(0) << "\"" << endl;
    ofs << "\"Restriction Enzyme: " << _cmd.GetEndonuclease(0) << "\"" << endl;
    ofs << "\"Query allowed at most " << _cmd.Mismatch() << " mismatches within ";
    ofs << _cmd.MaxBase() << " bases from 5\' end of primer.\"" << endl << endl;
    ofs << "\"Forward\",\"Reverse\",\"Forward\",\"Reverse\",\"Abundance\",\"Name\"" << endl;

    for (list<stNICHE>::iterator i = _niche.begin(); !(i == _niche.end()); ++i)
    {
        sprintf(buffer, "%.0f,%.0f,%.2f,%.2f,%.6f,\"%s\"",
            (*i).fpredict,            // predicted forward fragment
            (*i).rpredict,            // predicted reverse fragment
            (*i).fobserve,            // observed (sample) forward fragment
            (*i).robserve,            // observed (sample) reverse fragment
            (*i).biomass,             // normalized relative abundance
            (*i).organism.c_str());  // name of the species
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteCSV()

/*
 * draw the output in PHP format for web display
*/
bool WritePHP(
    CmdParam& _cmd)    // command-line parameters
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
    ofs << "  DrawHeader(\"MiCA: T-RFLP Analysis (APLAUS+) Output\");" << endl;
    ofs << "  DrawTRFLP(" << _cmd.GetFilename() << ");" << endl;
    ofs << "  DrawFooter();" << endl << "?>" << endl;
    ofs.close(); return(true);
}   // end of WritePHP()

// globally accessible classes for multithreading
CmdParam cmd; SeqDB rdp;
list<stNICHE> niche;
pthread_mutex_t mtxLock;    // critical region lock for database

/*
 * the prodcution trflp function
*/
void* DoTRFLP(void*)
{
    int forward, reverse; stNICHE item; bool run = false;

    pthread_mutex_lock(&mtxLock);
    // ** enter the critical section for records
    tRFLP rflp(cmd);
    // ** leave the critical section for database
    pthread_mutex_unlock(&mtxLock);

    do
    {
        pthread_mutex_lock(&mtxLock);
        // ** enter the critical section for database
        run = rdp.NextRecord();     // retrive a sequence from the database

        if (run)
        {
            rflp.SetStrand(rdp.GetOrigin());
            item.organism = rdp.GetOrganism();
        }   // set the sequence for the search
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLock);

        if (!run || !rflp.Delimit())
        {
            continue;
        }   // skip if there is no more sequences or amplification fails

        rflp.Digest(forward, reverse);    // perform restriction digest
        item.fpredict = static_cast<double>(forward),
        item.rpredict = static_cast<double>(reverse);

        if (!rflp.MatchSample(item))
        {   
            continue;
        }   // both fragments must match to be included in the list

        pthread_mutex_lock(&mtxLock);
        // ** enter the critical section for records
        niche.push_back(item);
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLock);
    } while (run);

    return(NULL);
}   // end of DoTRFLP()

/*
 * to compile, type:
 * g++ cmdparam.cpp seqdb.cpp bitvector.cpp trflp.cpp -o trflp.exe
*/
int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "insufficient number of parameters" << endl;
        return(1);
    }

    cmd.OpenFile(argv[1]);                  // open the parameter file
    rdp.OpenFile(cmd.GetDatabase());        // open the sequence database
    niche.clear();

    pthread_mutex_init(&mtxLock, NULL);     // initialize the lock for database
    vector<pthread_t> pts(nMaxTHREAD, 0);   // a vector for pthread

    for (unsigned int i = 0; i < pts.size(); ++i)
    {
        pthread_create(&pts[i], NULL, &DoTRFLP, NULL);
    }   // first, spawn the threads

    for (unsigned int j = 0; j < pts.size(); ++j)
    {
        pthread_join(pts[j], NULL);
    }   // wait for all threads to complete

    bool (*SortOption[8])(const stNICHE&, const stNICHE&) =
    {
        SortForwardA,   // sort by the sample forward fragment size in ascending order
        SortReverseA,   // sort by the sample reverse fragment size in ascending order
        SortBiomassA,   // sort by the species abundnace in ascending order
        SortOrganismA,  // sort by the species name in ascending order
        SortForwardD,   // sort by the sample forward fragment size in descending order
        SortReverseD,   // sort by the sample reverse fragment size in descending order
        SortBiomassD,   // sort by the species abundance in descending order
        SortOrganismD   // sort by the species name is descending order
    };  // nasty function pointers

    tRFLP rflp(cmd);
    rflp.SetAbundance(niche);     // calculate the relative abundance
    niche.sort(SortOption[cmd.SortOption()]);

    /*
     * write the output in various formats; explicitly signal the compiler that these
     * funcation calls do not require any specific order. the compiler is free to
     * rearrange for processor dispatch and parallel processing.
    */
    WriteCSV(niche, cmd), WriteTXT(niche, cmd), WritePHP(cmd), WriteDAT(niche, cmd);

    pthread_exit(NULL);
    pthread_mutex_destroy(&mtxLock);
    return(0);
}   // end of main()
