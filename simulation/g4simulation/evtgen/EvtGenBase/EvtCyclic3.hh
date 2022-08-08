
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#ifndef EVT_CYCLIC3_HH
#define EVT_CYCLIC3_HH

#include <iosfwd>

// Cyclic permutations of three indices A,B,C and their parings

namespace EvtCyclic3 {

    enum Index
    {
        A = 0,
        B = 1,
        C = 2
    };
    enum Pair
    {
        BC = 0,
        CB = BC,
        CA = 1,
        AC = CA,
        AB = 2,
        BA = AB
    };
    enum Perm
    {
        ABC = 0,
        BCA = 1,
        CAB = 2,
        CBA = 3,
        BAC = 4,
        ACB = 5
    };

    // Permutations (multiplication is not transitive)

    Index permute( Index i, Perm p );
    Perm permutation( Index i1, Index i2, Index i3 );
    Perm permute( Perm i, Perm p );
    Pair permute( Pair i, Perm p );

    Pair i2pair( int i );

    // Index-to-index

    Index prev( Index i );
    Index next( Index i );
    Index other( Index i, Index j );

    // Index-to-pair

    Pair other( Index i );
    Pair combine( Index i, Index j );

    // Pair-to-pair conversions

    Pair prev( Pair i );
    Pair next( Pair i );
    Pair other( Pair i, Pair j );

    // Pair-to-index conversions

    Index first( Pair i );
    Index second( Pair i );
    Index other( Pair i );
    Index common( Pair i, Pair j );

    // String to Index, Pair

    Index strToIndex( const char* str );
    Pair strToPair( const char* str );

    // To string conversions

    const char* c_str( Index i );
    const char* c_str( Pair i );
    const char* c_str( Perm i );

    // Useful name strings

    char* append( const char* str, EvtCyclic3::Index i );
    char* append( const char* str, EvtCyclic3::Pair i );

}    // namespace EvtCyclic3

//where should these go?
//ostream& operator<<(ostream&, EvtCyclic3::Index);
//ostream& operator<<(ostream&, EvtCyclic3::Pair);

#endif
