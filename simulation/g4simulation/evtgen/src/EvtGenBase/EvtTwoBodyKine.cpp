
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

#include "EvtGenBase/EvtTwoBodyKine.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::endl;
using std::ostream;

EvtTwoBodyKine::EvtTwoBodyKine() : _mA( 0. ), _mB( 0. ), _mAB( 0. )
{
}

EvtTwoBodyKine::EvtTwoBodyKine( double mA, double mB, double mAB ) :
    _mA( mA ), _mB( mB ), _mAB( mAB )
{
    if ( mAB < mA + mB ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << mAB << " < " << mA << " + " << mB << endl;
        assert( 0 );
    }
}

double EvtTwoBodyKine::m( Index i ) const
{
    double ret = _mAB;
    if ( A == i )
        ret = _mA;
    else if ( B == i )
        ret = _mB;

    return ret;
}

double EvtTwoBodyKine::p( Index i ) const
{
    double p0 = 0.;

    if ( i == AB ) {
        double x = _mAB * _mAB - _mA * _mA - _mB * _mB;
        double y = 2 * _mA * _mB;
        p0 = sqrt( x * x - y * y ) / 2. / _mAB;
    } else if ( i == A ) {
        double x = _mA * _mA - _mAB * _mAB - _mB * _mB;
        double y = 2 * _mAB * _mB;
        p0 = sqrt( x * x - y * y ) / 2. / _mA;
    } else {
        double x = _mB * _mB - _mAB * _mAB - _mA * _mA;
        double y = 2 * _mAB * _mA;
        p0 = sqrt( x * x - y * y ) / 2. / _mB;
    }

    return p0;
}

double EvtTwoBodyKine::e( Index i, Index j ) const
{
    double ret = m( i );
    if ( i != j ) {
        double pD = p( j );
        ret = sqrt( ret * ret + pD * pD );
    }
    return ret;
}

void EvtTwoBodyKine::print( ostream& os ) const
{
    os << " mA = " << _mA << endl;
    os << " mB = " << _mB << endl;
    os << "mAB = " << _mAB << endl;
}

ostream& operator<<( ostream& os, const EvtTwoBodyKine& p )
{
    p.print( os );
    return os;
}
