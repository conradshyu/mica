/*
 * SEQDB.H
 *
 * Written by Conrad Shyu (shyu4751@uidaho.edu)
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
 * this program provides an uninform interface to access the sequence database
 * this file is required for all operations
 *
 * All rights reserved. Copyright (R) 2004.
 * last updated on June 26, 2004
 * revised on December 24, 2009
*/
#ifndef _SEQDB_H
#define _SEQDB_H

#include <list>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

/*
 * class implementation to parse the sequences from the plain text format
*/
class   SeqDB
{
public:
    SeqDB() {};     // default constructor
    SeqDB(const char*);
    ~SeqDB() { ifInFile.close(); }

    // inline functions
    const string& GetLocus() const      { return(szLocus); }
    const string& GetOrigin() const     { return(szOrigin); }
    const string& GetOrganism() const   { return(szOrganism); }
    const string& GetAccession() const  { return(szAccession); }

    bool OpenFile(const char*);
    bool NextRecord();          // retrieve the next available sequence
    void PrintLocus();
    void PrintOrigin();
    void PrintOrganism();
    void PrintAccession();

private:
    ifstream ifInFile;
    string szAccession, szLocus, szOrigin, szOrganism;
};

#endif  // _SEQDB_H
