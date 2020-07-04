/*
 * BITVECTOR.H
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
 * this program implements the encoding scheme for the sequences
 * this file is required for all operations
 *
 * All rights reserved. Copyright (R) 2005.
 * last updated on April 15, 2005
*/
#ifndef _BITVECTOR_H
#define _BITVECTOR_H

#include <bitset>
#include <cctype>
#include <string>
#include <iostream>

using namespace std;

const int nMaxINT_WIDTH     = 8 * sizeof(int);
const int nMaxNUCLEOTIDE    = 4;
enum { baseA = 0, baseC, baseG, baseT };

/*
 * class implementation to convert the character sequences into binary stream
 * currently the length of the binary stream is limited to 32 bits. since almost
 * all primers are less than 25 base pairs long, this restriction should not be
 * a problem
*/
class   BitVector
{
public:
    BitVector(const int = 0, const int = 0);
    BitVector(const BitVector& _bv)   { *this = _bv; }
    ~BitVector() {};

    BitVector& operator=(const BitVector&);
    unsigned int operator&(const BitVector&) const;   // AND two bitvectors
    unsigned int operator^(const BitVector&) const;   // XOR two bitvectors
    unsigned int operator[](const int) const;

    int SetForwardPrimer(const string&);
    int SetReversePrimer(const string&);
    int SetEndonuclease(const string&);
    int SetForwardStrand(int, const string&);
    int SetReverseStrand(int, const string&);
    int SetDigestStrand(int, const string&);
    void AddNucleotide(char);
    void SetMismatch(const int = 0, const int = 0);
    void GetOffset(int&, int&) const;

    int GetOffset() const               { return(nOffset); }
    int GetDistance() const             { return(szStrand.length() - nDistance); }
    int GetLength() const               { return(szStrand.length()); }
    unsigned int GetMask() const        { return(uMask); }
    unsigned int GetConserved() const   { return(uConserved); }
    const string& GetStrand() const     { return(szStrand); }

    unsigned int GetAdenine() const  { return((*this)[baseA]); }
    unsigned int GetCytosine() const { return((*this)[baseC]); }
    unsigned int GetGuanine() const  { return((*this)[baseG]); }
    unsigned int GetThymine() const  { return((*this)[baseT]); }

    bool IsEnzyme(const BitVector&);
    bool IsPrimer(const BitVector&);
    void Print();       // print out the current bit configurations
    void Clear();       // clear the bit patterns in the class

private:
    unsigned int uStrand[nMaxNUCLEOTIDE];
    unsigned int uMask, uConserved;
    string szStrand;
    int nLeftOffset, nRightOffset, nOffset;
    int nConserved, nExact, nMaxBase, nDistance;

    void SetContrast();         // calculate the bit patterns for mismatch
};  // end of class definition for BitVector

#endif  // _BITVECTOR_H
