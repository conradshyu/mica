/*
 * BITVECTOR.CPP
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
 * NOTE: although the methods SetForwardPrimer(), SetReversePrimer(),
 * SetForwardTemplate(), and SetReverseTemplate() look very similar, however there
 * are subtle and critical differences between them
 *
 * All rights reserved. Copyright (R) 2005.
 * last updated on April 15, 2005
 * last updated on July 7, 2007
*/
#include <bitvector.h>

// for debugging purpose
//#define _VERBOSE

/*
 * note: the definition of the following bit patterns dictates how the ambiguities are
 * handled in the search algorithm. the bit patterns significantly simplies and speeds
 * up the translations of characters to binary representations.
*/

// ambiguity patterns for the template sequences, N and X are not searched
const unsigned int uSource[nMaxNUCLEOTIDE] =
{
    //                             ZY XWVU TSRQ PONM LKJI HGFE DCBA
    0x00621089, // adenine  0000:0000:0110:0010:0001:0000:1000:1001
    0x01241086, // cytosine 0000:0001:0010:0100:0001:0000:1000:0110
    0x0026044a, // guanine  0000:0000:0010:0110:0000:0100:0100:1010
    0x0158048a  // thymine  0000:0001:0101:1000:0000:0100:1000:1010
};

/*
 * support of inosine (I), inosine basically translates to N, resolve to A, C, G, and T
 * added by Conrad Shyu, July 7, 2007
 *
 * Inosine is a nucleotide that is formed when hypoxanthine is attached to a ribose ring
 * (also known as a ribofuranose) via a beta-Ng-glycosidic bond. It is commonly found in
 * tRNAs and is essential for proper translation of the genetic code in wobble base
 * pairs. Clinical studies have shown that inosine has neuroprotective properties. It
 * has been proposed for administration after stroke, because observation has show that
 * axonal rewiring is encouraged. It has also been tried for multiple sclerosis and is
 * currently in phase II of the trials.
 *
 * The support of inosine was suggested by Roey Angel of Arava Insititute for
 * Environmental Studies on July 2, 2007
*/

// ambiguity patterns for the primer and endonuclease sequences, N and X are searched
const unsigned int uTarget[nMaxNUCLEOTIDE] =
{
    //                             ZY XWVU TSRQ PONM LKJI HGFE DCBA
    0x00E23189, // adenine  0000:0000:1110:0010:0011:0001:1000:1001
    0x01a43186, // cytosine 0000:0001:1010:0100:0011:0001:1000:0110
    0x00a6254a, // guanine  0000:0000:1010:0110:0010:0101:0100:1010
    0x01d8258a  // thymine  0000:0001:1101:1000:0010:0101:1000:1010
};

/*
 * default constructor; just initialize the variables
*/
BitVector::BitVector(
    const int _mis, const int _max)
{
    uMask = 0;
    nLeftOffset = nRightOffset = 0;
    SetMismatch(_mis, _max); Clear();
}   // end of class constructor

/*
 * overload the assignment operator; used as a copy constructor
*/
BitVector& BitVector::operator=(
    const BitVector& _bv)
{
    if (!(this == &_bv))
    {
        uStrand[baseA] = _bv.GetAdenine(); uStrand[baseC] = _bv.GetCytosine();
        uStrand[baseG] = _bv.GetGuanine(); uStrand[baseT] = _bv.GetThymine();
        uMask = _bv.GetMask(); szStrand = _bv.GetStrand();
        _bv.GetOffset(nLeftOffset, nRightOffset);
    }   // make sure don't copy itself

    return(*this);
}   // end of operator overload for =

/*
 * retrieve the bit streams for specified nucleotide
*/
unsigned int BitVector::operator[](
    const int _idx) const
{
    if (!(_idx < nMaxNUCLEOTIDE))
    {
        return(0);
    }   // make sure the index is within the range

    return(uStrand[_idx]);
}   // end of operator overload for []

/*
 * overload the 'AND' operator to work on four bitset template variables
*/
unsigned int BitVector::operator&(
    const BitVector& _bv) const
{
    unsigned int a = 0;

    a |= (uStrand[baseA] & _bv[baseA]);
    a |= (uStrand[baseC] & _bv[baseC]);
    a |= (uStrand[baseG] & _bv[baseG]);
    a |= (uStrand[baseT] & _bv[baseT]);

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(a);
    cout << "AND operation: " << bs << endl;
#endif  // _VERBOSE

    return(a);
}   // end of operator overload for &

/*
 * overload the 'XOR' operator to work on four bitset template variables
*/
unsigned int BitVector::operator^(
    const BitVector& _bv) const
{
    unsigned int a = 0;

    a |= (uStrand[baseA] ^ _bv[baseA]);
    a |= (uStrand[baseC] ^ _bv[baseC]);
    a |= (uStrand[baseG] ^ _bv[baseG]);
    a |= (uStrand[baseT] ^ _bv[baseT]);

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(a);
    cout << "XOR operation: " << bs << endl;
#endif  // _VERBOSE

    return(a);
}   // end of operator overload ^

/*
 * convert the forward primer into bit streams
 * note: this implementation gets rid of the long, ugly switch/case statement,
 * and improves the performance
*/
int BitVector::SetForwardPrimer(
    const string& _primer)
{
    szStrand = _primer; SetContrast();
    uMask = ~(0xFFFFFFFF << szStrand.length());   // mask to filter unwanted bits
    int base;

    for (unsigned int i = 0; i < szStrand.length(); ++i)
    {
        base = static_cast<int>(szStrand[i] - 'A');
        uStrand[baseA] <<= 0x1; uStrand[baseC] <<= 0x1;     // advance one bit
        uStrand[baseG] <<= 0x1; uStrand[baseT] <<= 0x1;

        // now, "add" a nucleotide into the matrix
        uStrand[baseA] |= ((uTarget[baseA] >> base) & 0x1);
        uStrand[baseC] |= ((uTarget[baseC] >> base) & 0x1);
        uStrand[baseG] |= ((uTarget[baseG] >> base) & 0x1);
        uStrand[baseT] |= ((uTarget[baseT] >> base) & 0x1);
    }   // convert characters into binary streams

    // mask out the unwanted bit streams
    uStrand[baseA] &= uMask; uStrand[baseC] &= uMask;
    uStrand[baseG] &= uMask; uStrand[baseT] &= uMask;

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(uMask);
    cout << "   uMask: " << bs << endl;
    cout << "szStrand: " << szStrand << endl;
#endif  // _VERBOSE

    return(szStrand.length());
}   // end of SetForwardPrimer()

/*
 * conver the reverse primer into bit streams
 * complement: A<->T, C<->G
*/
int BitVector::SetReversePrimer(
    const string& _primer)
{
    szStrand = _primer; SetContrast();
    uMask = ~(0xFFFFFFFF << szStrand.length());   // mask to filter unwanted bits
    int base;

    for (unsigned int i = 0; i < szStrand.length(); ++i)
    {
        base = static_cast<int>(szStrand[i] - 'A');
        uStrand[baseA] <<= 0x1; uStrand[baseC] <<= 0x1;     // advance one bit
        uStrand[baseG] <<= 0x1; uStrand[baseT] <<= 0x1;

        // now, "add" a nucleotide into the matrix
        // NOTE: watch out for the complement, A<->T, C<->G
        uStrand[baseA] |= ((uTarget[baseT] >> base) & 0x1);
        uStrand[baseC] |= ((uTarget[baseG] >> base) & 0x1);
        uStrand[baseG] |= ((uTarget[baseC] >> base) & 0x1);
        uStrand[baseT] |= ((uTarget[baseA] >> base) & 0x1);
    }   // convert characters into binary streams

    // mask out the unwanted bit streams
    uStrand[baseA] &= uMask; uStrand[baseC] &= uMask;
    uStrand[baseG] &= uMask; uStrand[baseT] &= uMask;

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(uMask);
    cout << "   uMask: " << bs << endl;
    cout << "szStrand: " << szStrand << endl;
#endif  // _VERBOSE

    return(szStrand.length());
}   // end of SetReversePrimer()

/*
 * convert the restriction enzyme into binary streams
 *
 * note: it has been noted that some users submitted their own restriction enzymes that
 * do not contain the cut '^' symbole. to resolve this issue, the cut site will be at the
 * first position of the restriction enzyme.
 * last updated on July 10, 2007
*/
int BitVector::SetEndonuclease(
    const string& _enzyme)
{
    szStrand.clear();
    uMask = ~(0xFFFFFFFF << (_enzyme.length() - 1));   // filter unwanted bits
    int base;
    nLeftOffset = 0;    // if not cut symbol '^' is found, the first position is assumed

    for (unsigned int i = 0; i < _enzyme.length(); ++i)
    {
        if (_enzyme[i] == '^')
        {
            nLeftOffset = i; continue;
        }   // skip the cut character

        szStrand += _enzyme[i];
        base = static_cast<int>(_enzyme[i] - 'A');
        uStrand[baseA] <<= 0x1; uStrand[baseC] <<= 0x1;     // advance one bit
        uStrand[baseG] <<= 0x1; uStrand[baseT] <<= 0x1;

        // now, "add" a nucleotide into the matrix
        uStrand[baseA] |= ((uTarget[baseA] >> base) & 0x1);
        uStrand[baseC] |= ((uTarget[baseC] >> base) & 0x1);
        uStrand[baseG] |= ((uTarget[baseG] >> base) & 0x1);
        uStrand[baseT] |= ((uTarget[baseT] >> base) & 0x1);
    }   // convert characters into binary streams

    // mask out the unwanted bit streams
    uStrand[baseA] &= uMask; uStrand[baseC] &= uMask;
    uStrand[baseG] &= uMask; uStrand[baseT] &= uMask;
    // calculate the cut offsets
    nRightOffset = szStrand.length() - nLeftOffset;

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(uMask);
    cout << "   uMask: " << bs << endl;
    cout << "szStrand: " << szStrand << endl;
#endif  // _VERBOSE

    return(szStrand.length());
}   // end of SetEndonuclease()

/*
 * convert the template sequence for the search of forward primer
 * note: N and X are not searched
*/
int BitVector::SetForwardStrand(
    int _length,                // length is dictated by the forward primer
    const string& _strand)     // template sequence
{
    szStrand.clear(); nDistance = _length;
    uMask = ~(0xFFFFFFFF << _length);     // mask to filter unwanted bits

    for (int i = 0; i < _length; ++i)
    {
        AddNucleotide(_strand[i]);
    }   // convert characters into binary streams

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(uMask);
    cout << "   uMask: " << bs << endl;
    cout << "szStrand: " << szStrand << endl;
#endif  // _VERBOSE

    return(szStrand.length());
}   // end of SetForwardStrand()

/*
 * convert the template sequence for the search of reverse primer
 * note: N and X are not searched
*/
int BitVector::SetReverseStrand(
    int _length,                // length is dictated by the reverse primer
    const string& _strand)     // template sequence
{
    szStrand.clear(); nDistance = _length;
    uMask = ~(0xFFFFFFFF << _length);     // mask to filter unwanted bits
    int back = _strand.length();

    for (int i = 0; i < _length; ++i)
    {
        AddNucleotide(_strand[--back]);
    }   // convert characters into binary streams

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(uMask);
    cout << "   uMask: " << bs << endl;
    cout << "szStrand: " << szStrand << endl;
#endif  // _VERBOSE

    return(szStrand.length());
}   // end of SetReverseStrand()

/*
 * convert the template sequence for the search of the restriction enzyme
*/
int BitVector::SetDigestStrand(
    int _length,                // length is dictated by the restriction enzyme
    const string& _strand)     // template sequence
{
    szStrand.clear(); nDistance = 0;
    uMask = ~(0xFFFFFFFF << _length);     // mask to filter unwanted bits

    for (int i = 0; i < _length; ++i)
    {
        AddNucleotide(_strand[i]);
    }   // convert characters into binary streams

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs(uMask);
    cout << "   uMask: " << bs << endl;
    cout << "szStrand: " << szStrand << endl;
#endif  // _VERBOSE

    return(szStrand.length());
}   // end of SetDigestStrand()

/*
 * add a nucleotide into the matrix
 * note: N and X are not searched
*/
void BitVector::AddNucleotide(
    char _base)                // the nucleotide to be added into the matrix
{
    szStrand += _base;

    int base = static_cast<int>(_base - 'A');
    uStrand[baseA] <<= 0x1; uStrand[baseC] <<= 0x1;     // advance one bit
    uStrand[baseG] <<= 0x1; uStrand[baseT] <<= 0x1;

    // now, "add" a nucleotide into the matrix
    // note: N and X are not searched
    uStrand[baseA] |= ((uSource[baseA] >> base) & 0x1);
    uStrand[baseC] |= ((uSource[baseC] >> base) & 0x1);
    uStrand[baseG] |= ((uSource[baseG] >> base) & 0x1);
    uStrand[baseT] |= ((uSource[baseT] >> base) & 0x1);

    // mask out the unwanted bit streams
    uStrand[baseA] &= uMask; uStrand[baseC] &= uMask;
    uStrand[baseG] &= uMask; uStrand[baseT] &= uMask;
}   // end of AddNucleotide()

/*
 * reset the number of mismatches allowed in the search
*/
void BitVector::SetMismatch(
    int _mis, int _max)
{
    nExact = _max - _mis; nMaxBase = _max;

#ifdef _VERBOSE
    cout << "  nExact: " << nExact << endl;
    cout << "nMaxBase: " << nMaxBase << endl;
#endif  // _VERBOSE
}   // end of SetMismatch()

/*
 * get the offset values
*/
void BitVector::GetOffset(
    int& _l, int& _r) const
{
    _l = nLeftOffset; _r = nRightOffset;
}   // end of GetOffset()

/*
 * check if there is a match on the restriction enzyme
 * it is assumed that the restriction enzymes are palindrome sequences
*/
bool BitVector::IsEnzyme(
    const BitVector& _bv)
{
    unsigned int found = 0;

    // search for the forward match on the restriction enzyme
    found |= uStrand[baseA] & _bv[baseA];
    found |= uStrand[baseC] & _bv[baseC];
    found |= uStrand[baseG] & _bv[baseG];
    found |= uStrand[baseT] & _bv[baseT];

#ifdef _VERBOSE
    bitset<32> bs;
    bs = found; cout << "enzyme match: " << bs << endl; bs.reset();
#endif  // _VERBOSE

    if (found == uMask)
    {
        nOffset = nRightOffset; return(true);
    }   // a match in the forward direction

    return(false);
}   // end of IsEnzyme()

/*
 * check if there a match on the primer
*/
bool BitVector::IsPrimer(
    const BitVector& _bv)
{
    unsigned int final = (*this & _bv);
    int count = 0;

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs;
    bs = final; cout << "final: " << bs << endl;
#endif  // _VERBOSE

    if ((~final & uConserved))
    {
        return(false);
    }   // first, check the conserved regions

    // now, check the mismatched regions
    final >>= nConserved;

    for (int i = 0; i < nMaxBase; ++i, final >>= 0x1)
    {
        count += (final & 0x1);
    }   // accumulate the number of matched bits

    if (count < nExact)
    {
        return(false);
    }   // mismatch is less then required

    return(true);     // otherwise, it should be a match
}   // end of IsPrimer()

/*
 * print out the bit streams for debugging purposes
 * note: I am too lazy to write a function to convert the bit streams into strings
 * for pretty printing. I will just the bitset template.
*/
void BitVector::Print()
{
    bitset<nMaxINT_WIDTH> bs;
    cout << "   " << szStrand << endl;
    bs.reset(); bs = uStrand[baseA]; cout << "A: " << bs << endl;
    bs.reset(); bs = uStrand[baseC]; cout << "C: " << bs << endl;
    bs.reset(); bs = uStrand[baseG]; cout << "G: " << bs << endl;
    bs.reset(); bs = uStrand[baseT]; cout << "T: " << bs << endl;
}   // end of Print(); debugging function

/*
 * set the bit patterns for fast detect of a match
*/
void BitVector::SetContrast()
{
    nConserved = szStrand.length() - nMaxBase;
    uConserved = (0x1 << nConserved) - 1;

#ifdef _VERBOSE
    bitset<nMaxINT_WIDTH> bs;
    bs = uConserved; cout << "nConserved: " << bs << endl;
#endif  // _VERBOSE
}   // end of SetContrast()

/*
 * clear the bit streams of the four nucleotides
*/
void BitVector::Clear()
{
    uStrand[baseA] = uStrand[baseC] = 0;
    uStrand[baseG] = uStrand[baseT] = 0;
}   // end of Clear()

/*
 * test driver program
*/
/*
int main()
{
    BitVector f(0, 0);
    BitVector t(0, 0);
    string a = "GAGTTTGATCRTGGCTCAGGATGAACGCTGKCGGTCTGCTYAACACATGCAAGTCGAACGAAAGT"
    "CTTCGGACTTAGTGACGGACGGGTGAGTAWCGCGTGAGAATCTGCCTCCAGGTCGGGGACAACAGTTGGAAACGGC"
    "TGCTAATCCCGGATAAGCCGAAAGGTAAAAGATTTATTGCCTGGAGAGGAGCTCGCGTCCGATTAGCTAGATGGTG"
    "AGGTAAGAGCTCACCATGGCGACGATCGGTAGCTGTTCTGAGAGGAAGATCAGCCACACTGGGACTGAGACACGGC"
    "CCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCGCAATGGGCGCAAGCCTGACGGAGCAAGACCGCGTGCGG"
    "GAGGAAGGCCCTTGGGTCGTAAACCGCTTTTCTTAGGGAAGAAGCTCTGACGGTACCTAAGGAATCAGCCTCGGCT"
    "AACTCCGTGCCAGCAGCCGCGGTAATACGGAGGAGGCAAGCGTTATCCGGAATCATTGGGCGTAAAGCGTCCGCAG"
    "GCGGCTAATCAAGTCTGCTGTCAAAGACTGGGGCTCAACCCTGGGAAGGCAGTGGAAACTGAGAAGCTAGAGTGCA"
    "GTAGGGGTAGAGGGAATTCCCAGTGTAGCGGTGAAATGCGTAGAGATTGGGAAGAACACCGGTGGCGAAAGCGCTC"
    "TACTGGGCTGTAACTGACACTGAGGGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCTAGCCG"
    "TAAACGATGGAAACTAGGCGTGGCTTGTATCGACCCGAGCCGTGCCGAAGCTAACGCGTTAAGTTTCCCGCCTGGG"
    "GAGTACGCACGCAAGTGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCG"
    "ATGCAACGCGAAGAACCTTACCAGGGCTTGACATCTGGCGAATCTTTCCGAAAAGAGAGAGTGCCTTAGGGAGCGC"
    "CAAGACAGGTGGTGCATGGCTGTCGTCAGCTCGYGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCT"
    "CGTCCTTAGTTGCCGCATTGAGTTGGGCACTCTAGGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACG"
    "TCAAGTCATCATGCCCCTTACGTCTTGGGCTACACACGTACTACAATGGTTGGGACAACGGGCAGCGAAGCGGCGA"
    "CGCCGAGCGAATCCCAGCAAACCCAGCCTCAGTTCAGATCGCAGGCTGCAACTCGCCTGCGTGAAGGAGGAATCGC"
    "TAGTAATCGCCGGTCAGCATACGGCGGTGAATCCGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGAAGT"
    "TGGCCACGCCCGAAGTCGCTACTCTAACCTTTCGAGGAGGAGGACGCCGAAGGCAGGGCTGATGACTGGGGTGAAG"
    "TCGTAACAAGGTAGCCGTACCGGAAGGTGTGGCTGGATCACC";
    t.SetReverseStrand(18, a); t.Print();
    f.SetReversePrimer("TCCCCTAGCTTTCGTCCC");
    int i = a.size() - 18; cout << "index: " << i << endl;

    for (int k = 0; k < 15; ++k)
    {
        t.AddNucleotide(a[--i]); t.Print();
    }

    for (int i = 0; i < 10; ++i)
    {
        t.Print();

        if (f.IsPrimer(t))
        {
            cout << "found" << endl; break;
        }

        t.AddNucleotide(a[40 - i]);
    }

}   // end of main()
*/
