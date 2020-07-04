/*
 * CMDPARAM.H
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
 * this is the header file for the parameter parser
 *
 * All rights reserved. Copyright (R) 2005.
 * last updated on April 15, 2005
*/
#ifndef _CMDPARAM_H
#define _CMDPARAM_H

#include <list>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

#define _PARAM_H

/*
 * class to handle the command-line parameters
*/
class   CmdParam
{
public:
    CmdParam() {};      // default constructor; do nothing
    CmdParam(const char*);
    ~CmdParam() { ifInFile.close(); }

    const string& GetForwardPrimer(int = 0);
    const string& GetReversePrimer(int = 0);
    const string& GetEndonuclease(int = 0);
    const list<string>& GetForwardPrimer()  { return(szForwardPrimer); }
    const list<string>& GetReversePrimer()  { return(szReversePrimer); }
    const list<string>& GetEndonuclease()   { return(szEndonuclease); }
    const char* GetFilename() const { return(szFilename.c_str()); }
    const char* GetDatabase() const { return(szDatabase.c_str()); }
    const char* GetForwardSample() const { return(szForwardSample.c_str()); }
    const char* GetReverseSample() const { return(szReverseSample.c_str()); }

    bool OpenFile(const char*);
    bool OutputAll() const          { return(bOutputAll); }
    bool OutputShort() const        { return(!bOutputAll); }
    int ForwardBin() const          { return(nForwardBin); }
    int ReverseBin() const          { return(nReverseBin); }
    int EndonucleaseCount() const   { return(szEndonuclease.size()); }
    int ForwardPrimerCount() const  { return(szForwardPrimer.size()); }
    int ReversePrimerCount() const  { return(szReversePrimer.size()); }
    int SortOption() const          { return(nSortOption); }
    int MaxBase() const             { return(nMaxBase); }
    int Mismatch() const            { return(nMismatch); }
    void Print();

private:
    list<string> szForwardPrimer, szReversePrimer, szEndonuclease;
    string szFilename, szDatabase;
    string szForwardSample, szReverseSample;
    ifstream ifInFile;

    int nSortOption, nMaxBase, nMismatch;
    int nForwardBin, nReverseBin;
    bool bOutputAll;

    bool Parse();       // parse the command-line parameter
};  // end of class definition for CmdParam

#endif  // _PARAM_H
