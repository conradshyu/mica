/*
 * PSPA.H
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
 * this program performs the digest with given forward and reverse primers, and
 * restriction enzymes
 * this program is required for all operation
 *
 * All rights reserved. Copyrights (R) 2005.
 * last updated on April 28, 2005
 * last revised on July 10, 2007
*/
#ifndef _PSPA_H
#define _PSPA_H

// C++ headers
#include <list>
#include <cstdio>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

// class implementations
#include "seqdb.h"
#include "cmdparam.h"
#include "bitvector.h"

// for debugging purpose
//#define _VERBOSE

using namespace std;

/*
 * base class for DNA manipulations and digestions
*/
class   cPSPA
{
public:
    cPSPA(CmdParam&);
    ~cPSPA() {};

    bool SetStrand(const string&);
    int Delimit(vector<bool>&, vector<bool>&);    // delimit sequences with two primers
    void PrintStrand() const    { cout << szStrand; }

private:
    // binary representations for the primers and restriction enzymes
    bool bForwardFound, bReverseFound;
    BitVector bvForwardStrand, bvReverseStrand;
    list<BitVector> bvForwardPrimer, bvReversePrimer;

    string szStrand;

    bool Search(list<BitVector>::iterator, list<BitVector>::iterator);
};  // end of class definition for cPSPA

/*
 * class constructor for primer sequence prevalence analysis
 * setup the bit patterns for forward and reverse primers
*/
cPSPA::cPSPA(
    CmdParam& _cmd)    // command-line parameters
{
    BitVector forward, reverse;
    bvForwardPrimer.assign(_cmd.ForwardPrimerCount(), forward);
    bvReversePrimer.assign(_cmd.ReversePrimerCount(), reverse);
    int fidx = 0; int ridx = 0;

    for (list<BitVector>::iterator f = bvForwardPrimer.begin();
        !(f == bvForwardPrimer.end()); ++f, ++fidx)
    {
        (*f).SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
        (*f).SetForwardPrimer(_cmd.GetForwardPrimer(fidx));

#ifdef _VERBOSE
        cout << "forward primer:" << endl; (*f).Print();
#endif  // _VERBOSE
    }   // convert the forward primers into bit streams

    for (list<BitVector>::iterator r = bvReversePrimer.begin();
        !(r == bvReversePrimer.end()); ++r, ++ridx)
    {
        (*r).SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
        (*r).SetReversePrimer(_cmd.GetReversePrimer(ridx));

#ifdef _VERBOSE
        cout << "reverse primer:" << endl; (*r).Print();
#endif  // _VERBOSE
    }
}   // end of class constructor

/*
 * set the sequence to be searched
*/
bool cPSPA::SetStrand(
    const string& _s)
{
    szStrand = _s;

    if (!(szStrand.length() > 0))
    {
        return(false);
    }

#ifdef _VERBOSE
    cout << "szStrand: " << szStrand << endl;
#endif

    return(true);
}   // end of SetStran()

/*
 * search for the forward and reverse primers
*/
bool cPSPA::Search(
    list<BitVector>::iterator   _fp,    // iterator of the forward primer
    list<BitVector>::iterator   _rp)   // iterator of the reverse primer
{
    int forward, reverse;

    // convert sequence into bit streams for forward primer search
    bvForwardStrand.SetForwardStrand((*_fp).GetLength(), szStrand);
    forward = (*_fp).GetLength() - 1;

#ifdef _VERBOSE
    cout << "forward strand:" << endl; bvForwardStrand.Print();
#endif  // _VERBOSE

    bvReverseStrand.SetReverseStrand((*_rp).GetLength(), szStrand);
    reverse = szStrand.length() - (*_rp).GetLength();

#ifdef _VERBOSE
    cout << "reverse strand:" << endl; bvReverseStrand.Print();
#endif  // _VERBOSE

    bForwardFound = bReverseFound = false;

    while (forward < reverse)
    {
        if (!bForwardFound)
        {
            bForwardFound = (*_fp).IsPrimer(bvForwardStrand);
            bvForwardStrand.AddNucleotide(szStrand[++forward]);
        }   // search for the forward primer

        if (!bReverseFound)
        {
            bReverseFound = (*_rp).IsPrimer(bvReverseStrand);
            bvReverseStrand.AddNucleotide(szStrand[--reverse]);
        }   // search for the reverse primer

        if (bForwardFound && bReverseFound)
        {
            return(true);
        }   // both primer must be found
    }   // search for the forward and reverse primer

    return(false);
}   // end of Search()

/*
 * restrict the sequences with two primers
*/
int cPSPA::Delimit(
    vector<bool>&   _f,     // indicate if the forward primer has been found
    vector<bool>&   _r)    // indicate if the reverse primer has been found
{
    for (list<BitVector>::iterator i = bvForwardPrimer.begin();
        !(i == bvForwardPrimer.end()); ++i)
    {
        for (list<BitVector>::iterator j = bvReversePrimer.begin();
            !(j == bvReversePrimer.end()); ++j)
        {
            Search(i, j);     // search for the forward and reverse primers
            _f.push_back(bForwardFound), _r.push_back(bReverseFound);
        }
    }   // iterate through the entire lists of forward and reverse primers

    return(_f.size());
}   // end of Delimit()

#endif  // _PSPA_H
