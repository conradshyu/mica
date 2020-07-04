/*
 * PAT.CPP
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
 * this is a modified version of PAT (phylogenetic assignment tool), initially developed
 * by Angela Kent; many people complained that PAT does not work properly and it is
 * cumbersome to generate a database from MiCA, and then import into PAT for analysis
 *
 * this program performs the terminal restriction fragment analysis (t-rflp) using
 * single fragment profile
 * this program is required for all operation
 *
 * All rights reserved. Copyrights (R) 2007.
*/

// support class implementation
#include <pat.h>
#include <pthread.h>

// force PHP to return immediately
#define CLOSE_PHP       { fclose(stdin); fclose(stdout); fclose(stderr); };

const unsigned int nMaxTHREAD = 4;

/*
 * sort by the species abundance in ascending order
*/
bool SortBiomassA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.biomass < _b.biomass));
}   // end of SortBiomassA()

/*
 * sort by the species abundance in descending order
*/
bool SortBiomassD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.biomass > _b.biomass));
}   // end of SortBiomassD

/*
 * sort by the species name in ascending order
*/
bool SortOrganismA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.organism < _b.organism));
}   // end of SortOrganismA()

/*
 * sort by the species name is descending order
*/
bool SortOrganismD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.organism > _b.organism));
}   // end of SortOrganismD()

/*
 * sort by the sample forward fragment size in ascending order
*/
bool SortFragmentA(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.predict < _b.predict));
}   // end of SortForwardA()

/*
 * sort by the sample forward fragment size in descending order
*/
bool SortFragmentD(
    const stNICHE& _a, const stNICHE& _b)
{
    return((_a.predict > _b.predict));
}   // end of SortForwardD()

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
        // write the fragment lengths into the file
        sprintf(buffer, "%.0f,%.2f,%.6f,\"%s\",\"%s\"",
            (*i).predict,             // predicted trflp forward fragments
            (*i).observe,             // observed (sample) forward fragments
            (*i).biomass,             // normalized relative abundance
            (*i).accession.c_str(),   // species genbank accession number
            (*i).organism.c_str());  // species name
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteDAT()

/*
 * write the output in the plain text format, only the shortest fragments
 * format:
 * predicted and observed forward fragments
 * normalized relative abundance
 * genbank accession number
 * organism name
 *
 * last updated on July 14, 2006
 * last updated on July 7, 2007
 * last updated on July 10, 2007
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

    ofs << "Query returned " << _niche.size() << " record(s)" << endl
        << "Forward Primer: " << _cmd.GetForwardPrimer(0) << ", "
        << "Reverse Primer: " << _cmd.GetReversePrimer(0) << endl
        << "Restriction Enzyme: " << _cmd.GetEndonuclease(0) << endl << endl
        << "Query allowed at most " << _cmd.Mismatch()
        << " mismatches within " << _cmd.MaxBase()
        << " bases from 5\' end of primer" << endl << endl;

    ofs << "Predicted Observed Abundance Accession Name" << endl;

    for (list<stNICHE>::iterator i = _niche.begin(); !(i == _niche.end()); ++i)
    {
        // write the fragment lengths into the file
        sprintf(buffer, "%9.0f %8.2f %9.6f %9s %s",
            (*i).predict,             // predicted trflp forward fragments
            (*i).observe,             // observed (sample) forward fragments
            (*i).biomass,             // normalized relative abundance
            (*i).accession.c_str(),   // species genbank accession number
            (*i).organism.c_str());  // species name
        ofs << buffer << endl;
    }   // iterate through the entire list and print out the contents

    ofs.close(); return(true);
}   // end of WriteTXT()

/*
 * write the output in CSV format; Excel and other database programs readable
 * format:
 * predicted forward and revesre fragments
 * normalized relative abundance
 * accession number
 * organism name
 *
 * last updated on July 14, 2006
 * last updated on July 7, 2007
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

    ofs << "\"Query returned " << _niche.size() << " record(s).\"" << endl
        << "\"Forward Primer: " << _cmd.GetForwardPrimer(0) << ", "
        << "Reverse Primer: " << _cmd.GetReversePrimer(0) << "\"" << endl
        << "\"Restriction Enzyme: " << _cmd.GetEndonuclease(0) << "\"" << endl
        << "\"Query allowed at most " << _cmd.Mismatch() << " mismatches within "
        << _cmd.MaxBase() << " bases from 5\' end of primer.\"" << endl << endl
        << "\"Predicted\",\"Observed\",\"Abundance\",\"Accession\",\"Name\"" << endl;

    for (list<stNICHE>::iterator i = _niche.begin(); !(i == _niche.end()); ++i)
    {
        // write the fragment lengths into the file
        sprintf(buffer, "%.0f,%.2f,%.6f,\"%s\",\"%s\"",
            (*i).predict,             // predicted trflp forward fragments
            (*i).observe,             // observed (sample) forward fragments
            (*i).biomass,             // normalized relative abundance
            (*i).accession.c_str(),   // species genbank accession number
            (*i).organism.c_str());  // species name
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

    ofs << "<?php" << endl
        << "  require \"shared.data.inc\";" << endl
        << "  DrawHeader(\"MiCA: Phylogenetic Analysis (PAT+) Output\");" << endl
        << "  DrawPAT(" << _cmd.GetFilename() << ");" << endl
        << "  DrawFooter();" << endl << "?>" << endl;
    ofs.close(); return(true);
}   // end of WritePHP()

// globally accessible classes for multithreading
CmdParam cmd; SeqDB rdp;
list<stNICHE> niche;
pthread_mutex_t mtxLock;    // critical region lock for database and records

/*
 * the prodcution trflp function
*/
void* DoTRFLP(void*)
{
    int forward, reverse; bool run = false; stNICHE item;

    pthread_mutex_lock(&mtxLock);
    // ** enter the critical section for database
    cPAT rflp(cmd);   // instantiate trflp class with parameters
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
            item.organism = rdp.GetOrganism(), item.accession = rdp.GetAccession();
        }   // set the sequence for the search
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLock);

        if (!run || !rflp.Delimit())
        {
            continue;           // if primers cannot be found, do nothing
        }   // delimit the sequences with two primers

        rflp.Digest(forward, reverse);    // perform restriction digest
        item.predict = static_cast<double>(forward);

        if (!rflp.MatchSample(item))
        {
            continue;
        }   // only use the forward fragment for species identification

        pthread_mutex_lock(&mtxLock);
        // ** enter the critical section for database
        niche.push_back(item);
        // ** leave the critical section for database
        pthread_mutex_unlock(&mtxLock);
    } while (run);

    return(NULL);
}   // end of DoTRFLP()

/*
 * to compile, type:
 * g++ cmdparam.cpp seqdb.cpp bitvector.cpp pat.cpp -o pat
 *
 * note: for better performance, type:
 * g++ -O3 -static cmdparam.cpp seqdb.cpp bitvector.cpp pat.cpp -o pat
 *
 * last updated on July 7, 2007
*/
int main(int argc, char** argv)
{
    cmd.OpenFile(argv[1]);              // open the parameter file
    rdp.OpenFile(cmd.GetDatabase());      // open the sequence database
    niche.clear();

    pthread_mutex_init(&mtxLock, NULL);   // initialize the lock for database
    vector<pthread_t> pts(nMaxTHREAD, 0);     // a vector for pthread

    for (unsigned int i = 0; i < pts.size(); ++i)
    {
        pthread_create(&pts[i], NULL, &DoTRFLP, NULL);
    }   // first, spawn the threads

    for (unsigned int j = 0; j < pts.size(); ++j)
    {
        pthread_join(pts[j], NULL);
    }   // wait for all threads to complete

    cPAT rflp(cmd);

    for (list<stNICHE>::iterator n = niche.begin(); !(n == niche.end()); n++)
    {
        rflp.MatchSample(*n);
    }   // match the sample again

    rflp.SetAbundance(niche);   // calculate the relative abundance

    bool (*SortOption[6])(const stNICHE&, const stNICHE&) =
    {
        SortFragmentA,  // sort by the sample forward fragment size in ascending order
        SortBiomassA,   // sort by the species abundnace in ascending order
        SortOrganismA,  // sort by the species name in ascending order
        SortFragmentD,  // sort by the sample forward fragment size in descending order
        SortBiomassD,   // sort by the species abundance in descending order
        SortOrganismD   // sort by the species name is descending order
    };  // nasty function pointers

    niche.sort(SortOption[cmd.SortOption()]);

    /*
     * write the output in various formats; explicitly signal the compiler that these
     * funcation calls do not require any specific order. the compiler is free to
     * rearrange for processor dispatch and parallel processing.
    */
    WriteTXT(niche, cmd), WritePHP(cmd), WriteCSV(niche, cmd), WriteDAT(niche, cmd);

    pthread_exit(NULL);
    pthread_mutex_destroy(&mtxLock);
    return(0);
}   // end of main()
