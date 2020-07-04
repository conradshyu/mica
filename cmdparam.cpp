/*
 * CMDPARAM.CPP
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
 * this program will open up the parameter file and parse it contents
 * this file is required for all operations
 *
 * All rights reserved. Copyrights 2004.
 * last updated on June 17, 2004
*/
#include <cmdparam.h>

// for debugging purpose
//#define _VERBOSE

const int nMaxBUFFER = 8192;
const char* szParamDELIMIT = " =\n\t";

/*
 * class constructor; get the parameter filename without extension
*/
CmdParam::CmdParam(
    const char* _szFile) : szFilename(_szFile)
{
    if (!OpenFile(_szFile))
    {
        exit(1);
    }   // make sure the parameter file can be opened and parsed
}   // class constructor

bool CmdParam::OpenFile(
    const char* _szFile)
{
    szFilename = _szFile; ifInFile.open(szFilename.c_str(), ios::in);

    // initialize some important variables
    nSortOption = nMaxBase = nMismatch = 0;
    nForwardBin = nReverseBin = 0;
    bOutputAll = true;

#ifdef _VERBOSE
    cout << "parameter filename: " << szInFile << endl;
#endif  // _VERBOSE

    if (!ifInFile)
    {
        cout << "fatal: parameter file cannot be opened" << endl;
        return(false);
    }   // open the parameter file

    if (!Parse())
    {
        cout << "fatal: failed to parse the parameter file" << endl;
        return(false);
    }   // make sure the file is in the right format

    return(true);
}   // end of OpenFile()

/*
 * get the forward primer with given index
*/
const string& CmdParam::GetForwardPrimer(
    int _idx)      // which forward primer?
{
    list<string>::iterator i;

    for (i = szForwardPrimer.begin(); !(i == szForwardPrimer.end()) && (_idx--); ++i)
        ;

#ifdef _VERBOSE
    cout << "forward primer: " << (*i) << endl;
#endif  // _VERBOSE

    return(*i);
}   // end of GetForwardPrimer()

/*
 * get the reverse primer with given index
*/
const string& CmdParam::GetReversePrimer(
    int _idx)      // which reverse primer?
{
    list<string>::iterator i;

    for (i = szReversePrimer.begin(); !(i == szReversePrimer.end()) && (_idx--); ++i)
        ;

#ifdef _VERBOSE
    cout << "reverse primer: " << (*i) << endl;
#endif  // _VERBOSE

    return(*i);
}   // end of GetReversePrimer();

/*
 * endonuclease: a restriction enzyme that cleaves a nucleic acid (DNA or RNA) at
 * specific internal sites in the nucleotide base sequence
 *
 * get the restriction enzyme with given index
*/
const string& CmdParam::GetEndonuclease(
    int _idx)      // which restriction enzyme?
{
    list<string>::iterator i;

    // iterate through the list with the iterator
    for (i = szEndonuclease.begin(); !(i == szEndonuclease.end()) && (_idx--); ++i)
        ;

#ifdef _VERBOSE
    cout << "restriction enzyme: " << (*i) << endl;
#endif  // _VERBOSE

    return(*i);
}   // end of GetEndonuclease()

/*
 * print the command-line parameters
*/
void CmdParam::Print()
{
    const char* option[10] =
    {
        "forward fragments in ascending order",
        "reverse fragments in ascending order",
        "shortest forward fragment in ascending order",
        "shortest reverse fragment in ascending order",
        "organism name in ascending order",
        "forward fragments in descending order",
        "reverse fragments in descending order",
        "shortest forward fragment in descending order",
        "shortest reverse fragment in descending order",
        "organism name in descending order"
    };

    // detailed descriptions on the parameter settings
    cout << "parameter & output filename: " << GetFilename() << endl;
    cout << "number of forward primer(s): " << ForwardPrimerCount() << endl;
    cout << "number of reverse primer(s): " << ReversePrimerCount() << endl;
    cout << "  number of endonuclease(s): " << EndonucleaseCount() << endl;
    cout << "                sort option: " << option[SortOption()] << endl;
    cout << "             using database: " << GetDatabase() << endl;
    cout << "at most " << Mismatch() << " mismatches is allowed within the first ";
    cout << MaxBase() << " base pairs" << endl;
    cout << "   forward sample fragments: " << GetForwardSample() << endl;
    cout << "   reverse sample fragments: " << GetReverseSample() << endl;
    cout << "       forward fragment bin: " << ForwardBin() << endl;
    cout << "       reverse fragment bin: " << ReverseBin() << endl;

    int i;

    for (i = 0; i < ForwardPrimerCount(); ++i)
    {
        cout << "forward primer (" << i << "): " << GetForwardPrimer(i) << endl;
    }

    for (i = 0; i < ReversePrimerCount(); ++i)
    {
        cout << "reverse primer (" << i << "): " << GetReversePrimer(i) << endl;
    }

    for (i = 0; i < EndonucleaseCount(); ++i)
    {
        cout << "endonuclease (" << i << "): " << GetEndonuclease(i) << endl;
    }
}   // end of Print(); debugging function

/*
 * parse the parameter file and extract search settings
*/
bool CmdParam::Parse()
{
    char buffer[nMaxBUFFER]; char* token;

    // keep reading the parameters file until it ends
    while (ifInFile.getline(buffer, nMaxBUFFER))
    {
        // ignore comments
        if (buffer[0] == '#')
        {
#ifdef _VERBOSE
            cout << "comment: " << buffer << endl;
#endif  // _VERBOSE
            continue;
        }

        token = strtok(buffer, szParamDELIMIT);

        // parse the keywords; additional keywords can be added here
        if (!(strcmp(token, "database")))
        {
            szDatabase = strtok(0, szParamDELIMIT);
        }
        else if (!(strcmp(token, "max_base")))
        {
            nMaxBase = atoi(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "mismatch")))
        {
            nMismatch = atoi(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "sort_option")))
        {
            nSortOption = atoi(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "forward_shift")))
        {
            nForwardBin = atoi(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "reverse_shift")))
        {
            nReverseBin = atoi(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "output_all")))
        {
            bOutputAll = static_cast< bool >(atoi(strtok(0, szParamDELIMIT)));
        }
        else if (!(strcmp(token, "forward")))
        {
            szForwardPrimer.push_back(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "reverse")))
        {
            szReversePrimer.push_back(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "enzyme")))
        {
            szEndonuclease.push_back(strtok(0, szParamDELIMIT));
        }
        else if (!(strcmp(token, "filename")))
        {
            szFilename = strtok(0, szParamDELIMIT);
        }
        else if (!(strcmp(token, "forward_sample")))
        {
            szForwardSample = strtok(0, szParamDELIMIT);
        }
        else if (!(strcmp(token, "reverse_sample")))
        {
            szReverseSample = strtok(0, szParamDELIMIT);
        }
        else
        {
#ifdef _VERBOSE
            cout << "unknown token: " << token << endl;
#endif  // _VERBOSE
            continue;       // igore any commands that are not recognized
        }
    }

    return(true);         // everything is fine
}   // end of Parse()

/*
 * test driver program; it assumes the parameter file named, test.param
*/
/*
int main(int argc, char** argv)
{
    CmdParam q; q.OpenFile(argv[1]); q.Print();
}   // end of main()
*/
