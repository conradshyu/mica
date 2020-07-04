/*
 * ERPA.H
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
 * this program performs the digest with given forward and reverse primers, and
 * restriction enzymes
 * this program is required for all operation
 *
 * All rights reserved. Copyrights (R) 2005.
 * last updated on May 3, 2005
 * last revised on July 10, 2007
*/
#ifndef _ERPA_H
#define _ERPA_H

// C++ headers
#include <list>
#include <cmath>
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
//#define _FRAGMENTS

using namespace std;

/*
 * base class for DNA manipulations and digestions
*/
class   cERPA
{
public:
    cERPA(CmdParam&);
    ~cERPA() {};

    bool SetStrand(const string&);
    bool Delimit();                 // delimit sequences with two primers
    bool Digest(const string&);   // cut sequences with restriction enzymes

    void GetFragment(int& _ff, int& _rf) const
    {
        _ff = nForwardFragment; _rf = nReverseFragment;
    }   // return the forward and reverse fragments

    void PrintStrand() const    { cout << szStrand; }

private:
    int nForwardIndex, nReverseIndex, nEnzymeIndex;
    int nForwardDistance, nReverseDistance;
    int nForwardFragment, nReverseFragment;
    bool bForwardFound, bReverseFound;

    // binary representations for the primers and restriction enzymes
    BitVector bvForwardPrimer, bvReversePrimer;
    BitVector bvForwardStrand, bvReverseStrand;
    BitVector bvEndonuclease, bvEnzymeStrand;

    string szStrand;
};  // end of class defintion for cERPA

/*
 * class constructor for trflp
 * setup the bit patterns for forward and reverse primers, and the restriction
 * enzymes
*/
cERPA::cERPA(
    CmdParam& _cmd)    // command-line parameters
{
    // convert the forward and reverse primers into bit streams
    bvForwardPrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvReversePrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvForwardPrimer.SetForwardPrimer(_cmd.GetForwardPrimer(0));
    bvReversePrimer.SetReversePrimer(_cmd.GetReversePrimer(0));

#ifdef _VERBOSE     // print out the forward and reverse primers
    cout << "forward primer: " << endl; bvForwardPrimer.Print();
    cout << "reverse primer: " << endl; bvReversePrimer.Print();
#endif  // _VERBOSE
}   // end of class constructor

/*
 * set the sequence to be searched
*/
bool cERPA::SetStrand(
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

    // convert primers and sequence into bit streams for search
    bvForwardStrand.SetForwardStrand(bvForwardPrimer.GetLength(), szStrand);
    bvReverseStrand.SetReverseStrand(bvReversePrimer.GetLength(), szStrand);
    nForwardIndex = bvForwardPrimer.GetLength() - 1;
    nReverseIndex = szStrand.length() - bvReversePrimer.GetLength();

#ifdef _VERBOSE
    cout << "forward strand: " << endl; bvForwardStrand.Print();
    cout << "reverse strand: " << endl; bvReverseStrand.Print();
#endif

    return(true);
}   // end of SetStrand()

/*
 * restrict the sequences with two primers
*/
bool cERPA::Delimit()
{
    bForwardFound = bReverseFound = false;
    nForwardDistance = nReverseDistance = 0;

    while (nForwardIndex < nReverseIndex)
    {
        if (!bForwardFound)
        {
            bForwardFound = bvForwardPrimer.IsPrimer(bvForwardStrand);
            bvForwardStrand.AddNucleotide(szStrand[++nForwardIndex]);
        }   // search for the forward primer

        if (!bReverseFound)
        {
            bReverseFound = bvReversePrimer.IsPrimer(bvReverseStrand);
            bvReverseStrand.AddNucleotide(szStrand[--nReverseIndex]);
        }   // search for the reverse primer

        if (bForwardFound && bReverseFound)
        {
            // pointer will one step further even though a match has been found
            nForwardDistance = bvForwardStrand.GetDistance() - 1;
            nReverseDistance = bvReverseStrand.GetDistance() - 1;
            szStrand = szStrand.substr(nForwardIndex, nReverseIndex - nForwardIndex + 1);

#ifdef _VERBOSE         // print out some crucial variables
            cout << "nForwardDistance: " << nForwardDistance << endl;
            cout << "nReverseDistance: " << nReverseDistance << endl;
            cout << "         szStand: " << szStrand << endl;
#endif  // _VERBOSE

            return(true);
        }   // if both the forward and reverse primers are found, record the distance
    }   // search for the forward and reverse primer

    return(false);        // both primers cannot be found
}   // end of Delimit()

/*
 * cut the sequence with the restriction enzyme(s)
*/
bool cERPA::Digest(
    const string& _enzyme)
{
    int prior, size, full;
    vector<int> trf;
    full = bvForwardPrimer.GetLength() + bvReversePrimer.GetLength() + szStrand.length();
    bvEndonuclease.SetEndonuclease(_enzyme);          // construct the bit patterns
    nEnzymeIndex = bvEndonuclease.GetLength() - 1;      // index starts from 0
    bvEnzymeStrand.SetDigestStrand(bvEndonuclease.GetLength(), szStrand);
    prior = 0;

    while (nEnzymeIndex < static_cast<int>(szStrand.length()))
    {
        if (bvEndonuclease.IsEnzyme(bvEnzymeStrand))    // found enzyme
        {
            size = bvEnzymeStrand.GetDistance() - bvEndonuclease.GetOffset() - prior;
            trf.push_back(size); prior += size; bvEnzymeStrand.Clear();

#ifdef _VERBOSE
            cout << "prior: " << prior << endl << " size: " << size << endl;
#endif  // _VERBOSE
        }

        bvEnzymeStrand.AddNucleotide(szStrand[++nEnzymeIndex]);
    }   // now, search for the restriction enzymes

    if (trf.empty())      // full length
    {
        nForwardFragment = nReverseFragment = full;
    }   // now, fix the first and last fragments
    else
    {
        nForwardFragment = trf.front() + bvForwardPrimer.GetLength();
        nReverseFragment = szStrand.length() - prior + bvReversePrimer.GetLength();
    }   // calculate the forward and reverse fragment sizes

#ifdef _VERBOSE
    cout << "forward length: " << nForwardFragment << endl;
    cout << "reverse length: " << nReverseFragment << endl;
#endif  // _VERBOSE

    // return true if the enzyme can be found; otherwise returns false
    return(!trf.empty());
}   // end of Digest()

#endif  // _ERPA_H
