/*
 * ISPAR.H
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
 * last updated on April 26, 2005
*/
#ifndef _ISPAR_H
#define _ISPAR_H

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
class   tRFLP
{
public:
    tRFLP(CmdParam&);
    ~tRFLP() {};

    bool SetStrand(const string&);
    bool Delimit();     // delimit sequences with two primers
    int Digest();       // cut sequences with restriction enzymes

    /*
     * forward distance: the distance between first nucleotide to forward primer
     * reverse distance: the distance between last nucleotide to reverse primer
    */
    void GetDistance(int& _fd, int& _rd) const
    {
        _fd = nForwardDistance, _rd = nReverseDistance;
    }

    /*
     * return the forward and reverse fragments
    */
    void GetFragment(vector<int>& _ff, vector<int>& _rf) const
    {
        _ff = vForwardFragment, _rf = vReverseFragment;
    }

    /*
     * return the shortest forward and reverse fragments
    */
    void GetFragment(int& _ff, int& _rf) const
    {
        _ff = nForwardShort, _rf = nReverseShort;
    }

    void PrintStrand() const    { cout << szStrand; }

private:
    int nForwardIndex, nReverseIndex, nEnzymeIndex;
    int nForwardShort, nReverseShort;
    int nForwardDistance, nReverseDistance;
    bool bForwardFound, bReverseFound;

    // binary representations for the primers and restriction enzymes
    BitVector bvForwardPrimer, bvReversePrimer;
    BitVector bvForwardStrand, bvReverseStrand, bvEnzymeStrand;
    list<BitVector> bvEndonuclease;
    vector<int> vForwardFragment, vReverseFragment;

    string szStrand;
};  // end of class definition for tRFLP

/*
 * class constructor for trflp
 * setup the bit patterns for forward and reverse primers, and the restriction
 * enzymes
*/
tRFLP::tRFLP(
    CmdParam& _cmd)    // command-line parameters
{
    BitVector enzyme;

    // convert the forward and reverse primers into bit streams
    bvForwardPrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvReversePrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvForwardPrimer.SetForwardPrimer(_cmd.GetForwardPrimer(0));
    bvReversePrimer.SetReversePrimer(_cmd.GetReversePrimer(0));

#ifdef _VERBOSE     // print out the forward and reverse primers
    cout << "forward primer: " << endl; bvForwardPrimer.Print();
    cout << "reverse primer: " << endl; bvReversePrimer.Print();
#endif  // _VERBOSE

    // convert the restriction enzymes into bit streams
    for (int i = 0; i < _cmd.EndonucleaseCount(); ++i)
    {
        enzyme.SetEndonuclease(_cmd.GetEndonuclease(i));
        bvEndonuclease.push_back(enzyme);

#ifdef _VERBOSE     // print out the restriction enzyme
        cout << "restriction enzyme: " << endl; enzyme.Print();
#endif  // _VERBOSE
    }
}   // end of class constructor

/*
 * set the sequence to be searched
*/
bool tRFLP::SetStrand(
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
bool tRFLP::Delimit()
{
    bForwardFound = bReverseFound = false;
    nForwardDistance = nReverseDistance = 0;

    // search for the forward and reverse primer
    while (nForwardIndex < nReverseIndex)
    {
        // search for the forward primer
        if (!bForwardFound)
        {
            bForwardFound = bvForwardPrimer.IsPrimer(bvForwardStrand);
            bvForwardStrand.AddNucleotide(szStrand[++nForwardIndex]);
        }

        // search for the reverse primer
        if (!bReverseFound)
        {
            bReverseFound = bvReversePrimer.IsPrimer(bvReverseStrand);
            bvReverseStrand.AddNucleotide(szStrand[--nReverseIndex]);
        }

        // if both the forward and reverse primers are found, record the distance
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
        }
    }

    return(false);        // both primers cannot be found
}   // end of Delimit()

/*
 * cut the sequence with the restriction enzyme(s)
*/
int tRFLP::Digest()
{
    int prior, size, full;
    vector<int> trf;
    vForwardFragment.clear(); vReverseFragment.clear();
    full = bvForwardPrimer.GetLength() + bvReversePrimer.GetLength() + szStrand.length();
    nForwardShort = nReverseShort = full;

    // loop through the list of restriction enzymes
    for (list<BitVector>::iterator i = bvEndonuclease.begin();
        !(i == bvEndonuclease.end()); ++i)
    {
        bvEnzymeStrand.SetDigestStrand((*i).GetLength(), szStrand);
        nEnzymeIndex = (*i).GetLength() - 1;
        prior = 0; trf.clear();

        // now, search for the restriction enzymes
        while (nEnzymeIndex < static_cast<int>(szStrand.length()))
        {
            if ((*i).IsEnzyme(bvEnzymeStrand))    // found enzyme
            {
                size = bvEnzymeStrand.GetDistance() - (*i).GetOffset() - prior;
                trf.push_back(size); prior += size; bvEnzymeStrand.Clear();

#ifdef _VERBOSE
                cout << "prior: " << prior << endl;
                cout << " size: " << size << endl;
#endif  // _VERBOSE
            }

            bvEnzymeStrand.AddNucleotide(szStrand[++nEnzymeIndex]);
        }

        // now, fix the first and last fragments
        if (trf.empty())      // full length
        {
            vForwardFragment.push_back(full);
            vReverseFragment.push_back(full);

#ifdef _VERBOSE
            cout << "full length: " << full << endl;
#endif  // _VERBOSE
        }
        else
        {
            size = trf.front() + bvForwardPrimer.GetLength();
            nForwardShort = (nForwardShort > size) ? size : nForwardShort;
            vForwardFragment.push_back(size);

#ifdef _VERBOSE
            cout << "forward length: " << size << endl;
#endif  // _VERBOSE

            size = szStrand.length() - prior + bvReversePrimer.GetLength();
            nReverseShort = (nReverseShort > size) ? size : nReverseShort;
            vReverseFragment.push_back(size);

#ifdef _VERBOSE
            cout << "reverse length: " << size << endl;
#endif  // _VERBOSE
        }
    }

    return(vForwardFragment.size());      // number of fragments in the list
}   // end of Digest()

#endif  // _TRFLP_H
