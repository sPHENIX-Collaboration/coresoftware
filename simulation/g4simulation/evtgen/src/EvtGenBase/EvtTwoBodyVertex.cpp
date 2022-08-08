
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

#include "EvtGenBase/EvtTwoBodyVertex.hh"

#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::endl;
using std::ostream;

// Default ctor can sometimes be useful

EvtTwoBodyVertex::EvtTwoBodyVertex() : _LL( 0 ), _p0( 0 )
{
}

EvtTwoBodyVertex::EvtTwoBodyVertex( double mA, double mB, double mAB, int L ) :
    _kine(), _LL( L ), _p0( 0 )
{
    // Kinematics is initialized only if the decay is above threshold

    if ( mAB > mA + mB ) {
        _kine = EvtTwoBodyKine( mA, mB, mAB );
        _p0 = _kine.p();
    }
}

EvtTwoBodyVertex::EvtTwoBodyVertex( const EvtTwoBodyVertex& other ) :
    _kine( other._kine ),
    _LL( other._LL ),
    _p0( other._p0 ),
    _f( ( other._f ) ? new EvtBlattWeisskopf( *other._f ) : nullptr )
{
}

EvtTwoBodyVertex& EvtTwoBodyVertex::operator=( const EvtTwoBodyVertex& other )
{
    _kine = other._kine;
    _LL = other._LL;
    _p0 = other._p0;
    _f.reset( other._f ? new EvtBlattWeisskopf( *other._f ) : nullptr );
    return *this;
}

void EvtTwoBodyVertex::set_f( double R )
{
    _f = std::make_unique<EvtBlattWeisskopf>( _LL, R, _p0 );
}

double EvtTwoBodyVertex::widthFactor( EvtTwoBodyKine x ) const
{
    assert( _p0 > 0. );

    double p1 = x.p();
    double ff = formFactor( x );
    double factor = pow( p1 / _p0, 2 * _LL + 1 ) * mAB() / x.mAB() * ff * ff;

    return factor;
}

double EvtTwoBodyVertex::phaseSpaceFactor( EvtTwoBodyKine x,
                                           EvtTwoBodyKine::Index i ) const
{
    double p1 = x.p( i );
    double factor = pow( p1, _LL );
    return factor;
}

double EvtTwoBodyVertex::formFactor( EvtTwoBodyKine x ) const
{
    double ff = 1.;

    if ( _f ) {
        double p1 = x.p();
        ff = ( *_f )( p1 );
    }

    return ff;
}

void EvtTwoBodyVertex::print( ostream& os ) const
{
    os << " mA = " << mA() << endl;
    os << " mB = " << mB() << endl;
    os << "mAB = " << mAB() << endl;
    os << "  L = " << _LL << endl;
    os << " p0 = " << _p0 << endl;
}

ostream& operator<<( ostream& os, const EvtTwoBodyVertex& v )
{
    v.print( os );
    return os;
}
