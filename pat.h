/*
 * PAT.H
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
 * this program performs the terminal restriction fragment analysis (t-rflp) using
 * single fragment profile
 * this program is required for all operation
 *
 * All rights reserved. Copyrights (R) 2007.
*/
#ifndef _PAT_H
#define _PAT_H

// C++ headers
#include <list>
#include <cmath>
#include <cstdio>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

// class implementations
#include <seqdb.h>
#include <cmdparam.h>
#include <bitvector.h>

// for debugging purpose
//#define _VERBOSE

using namespace std;

const int nMaxBUFFER = 8192;
const char* szValueDELIMIT = " ,;|\n\t";

// define the structure for the sample data
typedef struct
{
    double fragment;    // fragment size
    double biomass;     // abundance associate with the fragment
    unsigned int count; // number of fragments that has been matched
} stSAMPLE;

/*
 * this structure establishes a relational model to associate the sample fragments
 * with the predicted fragments
*/
typedef struct
{
    string organism;                // name of the organism that matched both fragments
    string accession;               // genbank accession number
    double observe;                 // matched sample forward terminal fragment
    double predict;                 // matched, predicted reverse fragment
    double biomass;                 // normalized, relative abundance
    list<stSAMPLE>::iterator index; // index to the item in the sample profile
} stNICHE;

/*
 * base class for DNA manipulations and digestions
*/
class   cPAT
{
public:
    cPAT(CmdParam&);
    ~cPAT() {};

    bool SetStrand(const string&);
    bool Delimit();                         // delimit sequences with two primers
    bool Digest(int&, int&);              // cut sequences with restriction enzymes
    bool MatchSample(stNICHE&);           // match the predicted and observed fragments
    bool SetAbundance(list<stNICHE>&);    // calculate the relative abundance of species

    void PrintStrand() const    { cout << szStrand; }

private:
    int nForwardIndex, nReverseIndex, nEnzymeIndex;
    bool bForwardFound, bReverseFound;
    double dForwardBin;
    list<stSAMPLE> lsSample;

    // binary representations for the primers and restriction enzymes
    BitVector bvForwardPrimer, bvReversePrimer, bvEndonuclease;
    BitVector bvForwardStrand, bvReverseStrand, bvEnzymeStrand;

    string szStrand;

    bool LoadSample(list<stSAMPLE>&, const char*);
};  // end of class definition tRFLP

/*
 * class constructor for trflp
 * setup the bit patterns for forward and reverse primers, and the restriction
 * enzymes
*/
cPAT::cPAT(
    CmdParam& _cmd)    // command-line parameters
{
    // convert the forward and reverse primers into bit streams
    bvForwardPrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvReversePrimer.SetMismatch(_cmd.Mismatch(), _cmd.MaxBase());
    bvForwardPrimer.SetForwardPrimer(_cmd.GetForwardPrimer(0));
    bvReversePrimer.SetReversePrimer(_cmd.GetReversePrimer(0));
    bvEndonuclease.SetEndonuclease(_cmd.GetEndonuclease(0));
    dForwardBin = static_cast<double>(_cmd.ForwardBin());

#ifdef _VERBOSE     // print out the forward and reverse primers
    cout << "     forward primer: " << endl; bvForwardPrimer.Print();
    cout << "     reverse primer: " << endl; bvReversePrimer.Print();
    cout << " restriction enzyme: " << endl; bvEndonuclease.Print();
    cout << "forward window size: " << dForwardBin << endl;
#endif  // _VERBOSE

    // load the observed fragments into memory
    LoadSample(lsSample, _cmd.GetForwardSample());
}   // end of contructor

/*
 * set the sequence to be searched
*/
bool cPAT::SetStrand(
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
bool cPAT::Delimit()
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
bool cPAT::Digest(
    int&    _forward,       // predicted forward fragment
    int&    _reverse)      // predicted reverse fragment
{
    int size; int prior = 0;
    int full = bvForwardPrimer.GetLength() + bvReversePrimer.GetLength() + szStrand.length();
    bvEnzymeStrand.SetDigestStrand(bvEndonuclease.GetLength(), szStrand);
    nEnzymeIndex = bvEndonuclease.GetLength() - 1;      // index starts from 0
    vector<int> trf; trf.clear();                       // empty the fragment list

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
    }   // now, search for the restriction enzymes

    if (trf.empty())      // full length
    {
        _forward = _reverse = full;
    }   // now, fix the first and last fragments
    else
    {
        _forward = trf.front() + bvForwardPrimer.GetLength();
        _reverse = szStrand.length() - prior + bvReversePrimer.GetLength();
    }   // only the forward fragment is needed; reverse is not used

#ifdef _VERBOSE
    cout << "forward: " << _forward << ", reverse: " << _reverse << endl;
#endif  // _VERBOSE

    return(true);         // everything has worked as expected
}   // end of Digest()

/*
 * match the predicted fragments to the observed fragments in the sample file
*/
bool cPAT::MatchSample(
    stNICHE&    _niche)    // the structure for the plausible species
{
    for (list<stSAMPLE>::iterator i = lsSample.begin(); !(i == lsSample.end()); ++i)
    {
        if (!(fabs(_niche.predict - (*i).fragment) > dForwardBin))
        {
            _niche.observe = (*i).fragment, _niche.index = i, (*i).count++;
            return(true);
        }   // fragments are within the boundary
    }   // match the predicted and observed fragments

    return(false);
}   // end of MatchSample()

/*
 * load the user trflp profile into memory; the profile must be in the comma
 * separated value (csv) format
*/
bool cPAT::LoadSample(
    list<stSAMPLE>& _sample,    // structure of the fragments in the sample
    const char*     _name)     // name of the sample file
{
    ifstream ifs(_name, ios::in);
    char buffer[nMaxBUFFER]; stSAMPLE sample;

    if (!ifs)
    {
        return(false);
    }   // make sure the file can be opened successfully

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
 * normalize the relative abundance
*/
bool cPAT::SetAbundance(
    list<stNICHE>&  _niche)    // plausible community structure based on trflp data
{
    if (lsSample.empty() || _niche.empty())
    {
        return(false);
    }   // make sure both lists are not empty

    double abundance = 0.0;

    for (list<stSAMPLE>::iterator i = lsSample.begin(); !(i == lsSample.end()); ++i)
    {
        (*i).biomass /= ((*i).count == 0) ? 1.0 : static_cast<double>((*i).count);
    }   // first, normalize the abundance in the sample trflp profile

    for (list<stNICHE>::iterator j = _niche.begin(); !(j == _niche.end()); ++j)
    {
        (*j).biomass = (*(*j).index).biomass; abundance += (*j).biomass;
    }   // now, assign the abundance to the predicted trflp profile

    for (list<stNICHE>::iterator k = _niche.begin(); !(k == _niche.end()); ++k)
    {
        (*k).biomass /= abundance;

#ifdef _VERBOSE
        // print out the community profile
        cout << (*k).predict << ", " << (*k).observe << ", ";
        cout << (*k).biomass << ", " << (*k).accession << endl;
        cout << "name: " << (*k).organism << endl;
#endif  // _VERBOSE
    }   // finally, normalize the relative abundance for each species

    return(true);
}   // end of SetAbundance()

#endif  // _TRFLP_H
