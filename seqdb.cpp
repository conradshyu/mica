/*
 * SEQDB.CPP
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
 * All rights reserved. Copyright (R) 2005.
 * last updated on April 15, 2005
 * revised on December 24, 2009
*/
#include <seqdb.h>

// for debugging purpose
//#define _VERBOSE

const int   nDataBUFFER = 16384;
const char* szDataDELIMIT = "|\n";

/*
 * class constructor; open the sequence file
*/
SeqDB::SeqDB(
    const char* _szFile)
{
    if (!OpenFile(_szFile))
    {
        exit(1);
    }   // make sure the database file can be opened
}   // class constructor

bool SeqDB::OpenFile(
    const char* _szFile)
{
    ifInFile.open(_szFile, ios::in);

#ifdef _VERBOSE
    cout << "sequence filename: " << _szFile << endl;
#endif  // _VERBOSE

    if (!ifInFile)
    {
        cout << "cannot open sequence file: " << _szFile << endl;
        return(false);
    }   // make sure the database file can be opened

    return(true);
}   // end of OpenFile()

/*
 * look for the next record
*/
bool SeqDB::NextRecord()
{
    char buffer[nDataBUFFER];

    if (!ifInFile.getline(buffer, nDataBUFFER))
    {
        return(false);
    }   // there is no more elements in the cache buffer

    // now parse the contents; note: the order is critical
    szOrganism  = strtok(buffer, szDataDELIMIT);
    szAccession = strtok(0, szDataDELIMIT);
    szLocus     = strtok(0, szDataDELIMIT);
    szOrigin    = strtok(0, szDataDELIMIT);

#ifdef _VERBOSE
    cout << " szOrganism: " << szOrganism << endl;
    cout << "szAccession: " << szAccession << endl;
    cout << "    szLocus: " << szLocus << endl;
    cout << "   szOrigin: " << szOrigin << endl;
#endif  // _VERBOSE

    return(true);
}   // end of NextRecord()

/*
 * print out the entire list of locus names in the database
*/
void SeqDB::PrintLocus()
{
    while(NextRecord())
    {
        cout << szLocus << endl;
    }   // print out sequence locus information
}   // end of PrintLocus()

/*
 * print out the entire list of sequences in the database
 * note: the sequences can be extremely long
*/
void SeqDB::PrintOrigin()
{
    while(NextRecord())
    {
        cout << szOrigin << endl;
    }   // print out the sequence
}   // end of PrintOrigin()

/*
 * print out the entire list of organisms' names in the database
*/
void SeqDB::PrintOrganism()
{
    while(NextRecord())
    {
        cout << szOrganism << endl;
    }   // print out the name of organism
}   // end of PrintOrganism()

/*
 * print out the entire list of accession numbers in the database
*/
void SeqDB::PrintAccession()
{
    while(NextRecord())
    {
        cout << szAccession << endl;
    }   // print out access number
}   // end of PrintAccession()

/*
 * test driver program
*/
/*
int main(int argc, char** argv)
{
    SeqDB rdp; rdp.OpenFile("bacteria_all.txt");

    double mean = 0.0;
    double count = 0.0;

    while (rdp.NextRecord())
    {
        count += 1.0;
        mean += (static_cast<double>((rdp.GetOrigin()).size()) - mean) / count;
    }   // calculate the running average

    printf("total: %.0f, average: %.4f\n", count, mean);
}   // end of main()
*/
