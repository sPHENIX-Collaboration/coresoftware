
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

#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::ostream;

EvtVector4C::EvtVector4C()
{
    v[0] = EvtComplex( 0.0 );
    v[1] = EvtComplex( 0.0 );
    v[2] = EvtComplex( 0.0 );
    v[3] = EvtComplex( 0.0 );
}

EvtVector4C::EvtVector4C( const EvtComplex& e0, const EvtComplex& e1,
                          const EvtComplex& e2, const EvtComplex& e3 )
{
    v[0] = e0;
    v[1] = e1;
    v[2] = e2;
    v[3] = e3;
}

EvtVector4C rotateEuler( const EvtVector4C& rs, double alpha, double beta,
                         double gamma )
{
    EvtVector4C tmp( rs );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

EvtVector4C boostTo( const EvtVector4C& rs, const EvtVector4R p4 )
{
    EvtVector4C tmp( rs );
    tmp.applyBoostTo( p4 );
    return tmp;
}

EvtVector4C boostTo( const EvtVector4C& rs, const EvtVector3R boost )
{
    EvtVector4C tmp( rs );
    tmp.applyBoostTo( boost );
    return tmp;
}

void EvtVector4C::applyBoostTo( const EvtVector4R& p4 )
{
    double e = p4.get( 0 );

    EvtVector3R boost( p4.get( 1 ) / e, p4.get( 2 ) / e, p4.get( 3 ) / e );

    applyBoostTo( boost );

    return;
}

void EvtVector4C::applyBoostTo( const EvtVector3R& boost )
{
    double bx, by, bz, gamma, b2;

    bx = boost.get( 0 );
    by = boost.get( 1 );
    bz = boost.get( 2 );

    double bxx = bx * bx;
    double byy = by * by;
    double bzz = bz * bz;

    b2 = bxx + byy + bzz;

    if ( b2 == 0.0 ) {
        return;
    }

    assert( b2 < 1.0 );

    gamma = 1.0 / sqrt( 1 - b2 );

    double gb2 = ( gamma - 1.0 ) / b2;

    double gb2xy = gb2 * bx * by;
    double gb2xz = gb2 * bx * bz;
    double gb2yz = gb2 * by * bz;

    double gbx = gamma * bx;
    double gby = gamma * by;
    double gbz = gamma * bz;

    EvtComplex e2 = v[0];
    EvtComplex px2 = v[1];
    EvtComplex py2 = v[2];
    EvtComplex pz2 = v[3];

    v[0] = gamma * e2 + gbx * px2 + gby * py2 + gbz * pz2;

    v[1] = gbx * e2 + gb2 * bxx * px2 + px2 + gb2xy * py2 + gb2xz * pz2;

    v[2] = gby * e2 + gb2 * byy * py2 + py2 + gb2xy * px2 + gb2yz * pz2;

    v[3] = gbz * e2 + gb2 * bzz * pz2 + pz2 + gb2yz * py2 + gb2xz * px2;

    return;
}

void EvtVector4C::applyRotateEuler( double phi, double theta, double ksi )
{
    double sp = sin( phi );
    double st = sin( theta );
    double sk = sin( ksi );
    double cp = cos( phi );
    double ct = cos( theta );
    double ck = cos( ksi );

    EvtComplex x = ( ck * ct * cp - sk * sp ) * v[1] +
                   ( -sk * ct * cp - ck * sp ) * v[2] + st * cp * v[3];
    EvtComplex y = ( ck * ct * sp + sk * cp ) * v[1] +
                   ( -sk * ct * sp + ck * cp ) * v[2] + st * sp * v[3];
    EvtComplex z = -ck * st * v[1] + sk * st * v[2] + ct * v[3];

    v[1] = x;
    v[2] = y;
    v[3] = z;
}

ostream& operator<<( ostream& s, const EvtVector4C& v )
{
    s << "(" << v.v[0] << "," << v.v[1] << "," << v.v[2] << "," << v.v[3] << ")";

    return s;
}
