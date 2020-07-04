/*
 * TRFLP.H
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
 * this program performs the terminal restriction fragment analysis (t-rflp)
 * this program is required for all operation
 *
 * All rights reserved. Copyrights (R) 2005.
 * last updated on April 17, 2005
 * revised on November 15, 2005
 * major revision on the calculation of abundance on July 9, 2007
 * last revised on July 10, 2007
*/
#ifndef _TRFLP_H
#define _TRFLP_H

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

using namespace std;

const int nMaxBUFFER = 8192;
const char* szValueDELIMIT = " ,;|\n\t";

// define the structure for the sample data
typedef struct
{
    double fragment;        // fragment size
    double biomass;         // abundance associate with the fragment
    unsigned int count;     // number of fragments that has been matched
} stSAMPLE;

/*
 * this structure establishes a relational model to associate the sample fragments
 * with the predicted fragments
*/
typedef struct
{
    string organism;        // name of the organism that matched both fragments
    double fobserve;        // matched sample forward terminal fragment
    double robserve;        // matched sample reverse terminal fragment
    double fpredict;        // matched, predicted forward fragment
    double rpredict;        // matched, predicted reverse fragment
    double biomass;         // normalized, relative abundance
    list<stSAMPLE>::iterator findex;
    list<stSAMPLE>::iterator rindex;
} stNICHE;

/*
 * base class for DNA manipulations and digestions
*/
class   tRFLP
{
public:
    tRFLP(CmdParam&);
    ~tRFLP() {};

    bool SetStrand(const string&);
    bool Delimit();                         // delimit sequences with two primers
    bool Digest(int&, int&);              // cut sequences with restriction enzymes
    bool MatchSample(stNICHE&);           // match the predicted and observed fragments
    bool SetAbundance(list<stNICHE>&);    // calculate the relative abundance of species

    void PrintStrand() const    { cout << szStrand; }

private:
    int nForwardIndex, nReverseIndex, nEnzymeIndex;
    bool bForwardFound, bReverseFound;
    double dForwardBin, dReverseBin;

    // binary representations for the primers and restriction enzymes
    BitVector bvForwardPrimer, bvReversePrimer, bvEndonuclease;
    BitVector bvForwardStrand, bvReverseStrand, bvEnzymeStrand;
    list<stSAMPLE> lsForwardSample, lsReverseSample;

    string szStrand;

    bool LoadSample(list<stSAMPLE>&, const char*);
};  // end of class definition tRFLP

/*
 * class constructor for trflp
 * setup the bit patterns for forward and reverse primers, and the restriction
 * enzymes
*/
tRFLP::tRFLP(
    CmdParam& _cmd)    // command-line parameters
{
    // convert the forward and reverse primers into bit streams
    bvForwardPrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvReversePrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvForwardPrimer.SetForwardPrimer(_cmd.GetForwardPrimer(0));
    bvReversePrimer.SetReversePrimer(_cmd.GetReversePrimer(0));
    bvEndonuclease.SetEndonuclease(_cmd.GetEndonuclease(0));
    dForwardBin = static_cast<double>(_cmd.ForwardBin());
    dReverseBin = static_cast<double>(_cmd.ReverseBin());

#ifdef _VERBOSE     // print out the forward and reverse primers
    cout << "     forward primer: " << endl; bvForwardPrimer.Print();
    cout << "     reverse primer: " << endl; bvReversePrimer.Print();
    cout << " restriction enzyme: " << endl; bvEndonuclease.Print();
    cout << "forward window size: " << dForwardBin << endl;
    cout << "reverse window size: " << dReverseBin << endl;
#endif  // _VERBOSE

    // load the observed fragments into memory
    LoadSample(lsForwardSample, _cmd.GetForwardSample());
    LoadSample(lsReverseSample, _cmd.GetReverseSample());
}   // end of contructor

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
    }   // length of template sequence cannot be zero

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
            szStrand = szStrand.substr(nForwardIndex, nReverseIndex - nForwardIndex + 1);

#ifdef _VERBOSE         // print out some crucial variables
            cout << "szStand: " << szStrand << endl;
#endif  // _VERBOSE

            return(true);
        }   // if both the forward and reverse primers are found, record the distance
    }   // search for the forward and reverse primer

    return(false);        // both primers cannot be found
}   // end of Delimit()

/*
 * cut the sequence with the restriction enzyme(s)
*/
bool tRFLP::Digest(
    int&    _forward,   // predicted forward fragment
    int&    _reverse)  // predicted reverse fragment
{
    int size; int prior = 0;
    int full = bvForwardPrimer.GetLength() + bvReversePrimer.GetLength() + szStrand.length();
    bvEnzymeStrand.SetDigestStrand(bvEndonuclease.GetLength(), szStrand);
    nEnzymeIndex = bvEndonuclease.GetLength() - 1;      // index starts from 0
    vector<int> trf; trf.clear();       // empty the fragment list

    // now, search for the restriction enzymes
    while (nEnzymeIndex < static_cast<int>(szStrand.length()))
    {
        if (bvEndonuclease.IsEnzyme(bvEnzymeStrand))    // found enzyme
        {
            // calculate the digested fragment sizes
            size = bvEnzymeStrand.GetDistance() - bvEndonuclease.GetOffset() - prior;
            trf.push_back(size);          // save the digested fragment size
            bvEnzymeStrand.Clear();         // reset the bitvector for the next search
            prior += size;

#ifdef _VERBOSE
            cout << "prior: " << prior << endl;
            cout << " size: " << size << endl;
#endif  // _VERBOSE
        }

        bvEnzymeStrand.AddNucleotide(szStrand[++nEnzymeIndex]);
    }

    if (trf.empty())
    {
        _forward = _reverse = full;
    }   // full length, fix the first and last fragments
    else
    {
        _forward = trf.front() + bvForwardPrimer.GetLength();
        _reverse = szStrand.length() - prior + bvReversePrimer.GetLength();
    }   // calculate the forward and reverse fragment sizes

#ifdef _VERBOSE
    cout << "forward: " << _forward << ", reverse: " << _reverse << endl;
#endif  // _VERBOSE

    return(true);         // everything has worked as expected
}   // end of Digest()

/*
 * match the predicted fragments to the observed fragments in the sample file
*/
bool tRFLP::MatchSample(
    stNICHE&    _niche)    // predicted community profile
{
    if (lsForwardSample.empty() || lsReverseSample.empty())
    {
        return(false);
    }   // make sure both lists are not empty before matching proceeds

    bool is_forward = false;
    bool is_reverse = false;
    list<stSAMPLE>::iterator f, r;

    for (f = lsForwardSample.begin(); !(f == lsForwardSample.end()); ++f)
    {
        if (!(fabs((*f).fragment - _niche.fpredict) > dForwardBin))
        {
            _niche.fobserve = (*f).fragment, _niche.findex = f, is_forward = true;
            break;
        }   // forward fragment size is within the boundary
    }   // first, match the forward fragments

    for (r = lsReverseSample.begin(); !(r == lsReverseSample.end()); ++r)
    {
        if (!(fabs((*r).fragment - _niche.rpredict) > dReverseBin))
        {
            _niche.robserve = (*r).fragment, _niche.rindex = r, is_reverse = true;
            break;
        }   // reverse fragment size is within the boundary
    }   // then, match the reverse fragment

    if (is_forward && is_reverse)
    {
        (*f).count++, (*r).count++;
#ifdef _VERBOSE
        cout << "predicted and observed fragments are a match" << endl;
        cout << "(" << _niche.fpredict << ", " << _niche.fobserve << ") ";
        cout << "(" << _niche.rpredict << ", " << _niche.robserve << ")" << endl;
#endif  // _VERBOSE
    }   // if both fragments are a match, update the counter

    return(is_forward && is_reverse);
}   // end of MatchSample()

/*
 * load the user trflp profile into memory; the profile must be in the comma
 * separated value (csv) format
*/
bool tRFLP::LoadSample(
    list<stSAMPLE>& _sample,    // structure of the fragments in the sample
    const char* _name)         // name of the sample file
{
    ifstream ifs(_name, ios::in);
    char buffer[nMaxBUFFER]; stSAMPLE sample;

    if (!ifs)
    {
        return(false);
    }   // make sure the file can be opened or created successfully

    _sample.clear();        // empty the list before continue

    while (!ifs.eof())
    {
        ifs.getline(buffer, nMaxBUFFER);

        // each line should have at least three characters; "X,X"
        // kind of ugly, but works like a charm
        if (strlen(buffer) < 3)
        {
            continue;
        }

        // it is assumed that there are only two columns
        sample.biomass = atof(strtok(buffer, szValueDELIMIT));
        sample.fragment = atof(strtok(0, szValueDELIMIT));
        sample.count = 0;
        _sample.push_back(sample);        // save the record

#ifdef _VERBOSE
        // print out some information for debugging purpose
        cout << "fragment: " << sample.fragment << ", biomass: " << sample.biomass;
#endif  // _VERBOSE
    }   // read the file line by line until eof is encountered

    ifs.close(); return(true);
}   // end of LoadSample()

/*
 * calculate assign the normalized relative abundance to each species in the community
*/
bool tRFLP::SetAbundance(
    list<stNICHE>& _niche)
{
    if (_niche.empty())
    {
        return(false);
    }   // make sure the list is not empty

    list<stSAMPLE>::iterator i; list<stNICHE>::iterator j;
    double abundance = 0.0;

    for (i = lsForwardSample.begin(); !(i == lsForwardSample.end()); ++i)
    {
        (*i).biomass /= ((*i).count == 0) ? 1.0 : static_cast<double>((*i).count);
    }   // normalize the abundance for the forward fragments in the sample

    for (i = lsReverseSample.begin(); !(i == lsReverseSample.end()); ++i)
    {
        (*i).biomass /= ((*i).count == 0) ? 1.0 : static_cast<double>((*i).count);
    }   // normalize the abundance for the reverse fragments in the sample

    for (j = _niche.begin(); !(j == _niche.end()); ++j)
    {
        (*j).biomass = (*(*j).findex).biomass + (*(*j).rindex).biomass;
        abundance += (*j).biomass;        // accumulate the abundance
    }   // assign the normalized abundance to the predicted community profile

    for (j = _niche.begin(); !(j == _niche.end()); ++j)
    {
        (*j).biomass /= abundance;
    }   // finally, normalize the relative abundance in the predicted community profile

    return(true);
}   // end of SetAbundance()

#endif  // _TRFLP_H
