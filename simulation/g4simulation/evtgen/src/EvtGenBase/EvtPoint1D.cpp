
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

#include "EvtGenBase/EvtPoint1D.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <stdio.h>

EvtPoint1D::EvtPoint1D() :
    _min( 0. ), _max( -1. ), _value( 0. ), _valid( false )
{
}

EvtPoint1D::EvtPoint1D( double value ) :
    _min( 0. ), _max( -1. ), _value( value ), _valid( true )
{
}

EvtPoint1D::EvtPoint1D( double min, double max, double value ) :
    _min( min ),
    _max( max ),
    _value( value ),
    _valid( ( _min <= _value && _value <= _max ) ? true : false )
{
}

void EvtPoint1D::print() const
{
    printf( "%f (%f : %f)\n", _value, _min, _max );
}
