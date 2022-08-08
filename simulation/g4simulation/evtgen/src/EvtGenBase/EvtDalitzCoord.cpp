
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

#include "EvtGenBase/EvtDalitzCoord.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
using EvtCyclic3::Pair;
using std::endl;
using std::ostream;

// For coordinates it's good to alway have a
// default ctor. Initialize to something invalid.

EvtDalitzCoord::EvtDalitzCoord() :
    _i1( EvtCyclic3::AB ), _i2( EvtCyclic3::BC ), _q1( -1. ), _q2( -1. )
{
}

EvtDalitzCoord::EvtDalitzCoord( const EvtDalitzCoord& other ) :
    _i1( other._i1 ), _i2( other._i2 ), _q1( other._q1 ), _q2( other._q2 )
{
}

EvtDalitzCoord::EvtDalitzCoord( Pair i1, double q1, Pair i2, double q2 ) :
    _i1( i1 ), _i2( i2 ), _q1( q1 ), _q2( q2 )
{
}

bool EvtDalitzCoord::operator==( const EvtDalitzCoord& other ) const
{
    return ( _i1 == other._i1 && _i2 == other._i2 && _q1 == other._q1 &&
             _q2 == other._q2 );
}

void EvtDalitzCoord::print( ostream& os ) const
{
    os << _i1 << " " << _q1 << endl;
    os << _i2 << " " << _q2 << endl;
}

ostream& operator<<( ostream& os, const EvtDalitzCoord& p )
{
    p.print( os );
    return os;
}
