
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

#include "EvtGenBase/EvtVector3C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <iostream>
#include <math.h>
using std::ostream;

EvtVector3C::EvtVector3C()
{
    v[0] = EvtComplex( 0.0 );
    v[1] = EvtComplex( 0.0 );
    v[2] = EvtComplex( 0.0 );
}

EvtVector3C::EvtVector3C( const EvtComplex& e1, const EvtComplex& e2,
                          const EvtComplex& e3 )
{
    v[0] = e1;
    v[1] = e2;
    v[2] = e3;
}

EvtVector3C EvtVector3C::cross( const EvtVector3C& p2 )
{
    //Calcs the cross product.  Added by djl on July 27, 1995.

    EvtVector3C temp;

    temp.v[0] = v[1] * p2.v[2] - v[2] * p2.v[1];
    temp.v[1] = v[2] * p2.v[0] - v[0] * p2.v[2];
    temp.v[2] = v[0] * p2.v[1] - v[1] * p2.v[0];

    return temp;
}

EvtVector3C rotateEuler( const EvtVector3C& v, double alpha, double beta,
                         double gamma )
{
    EvtVector3C tmp( v );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

void EvtVector3C::applyRotateEuler( double phi, double theta, double ksi )
{
    EvtComplex temp[3];
    double sp, st, sk, cp, ct, ck;

    sp = sin( phi );
    st = sin( theta );
    sk = sin( ksi );
    cp = cos( phi );
    ct = cos( theta );
    ck = cos( ksi );

    temp[0] = ( ck * ct * cp - sk * sp ) * v[0] +
              ( -sk * ct * cp - ck * sp ) * v[1] + st * cp * v[2];
    temp[1] = ( ck * ct * sp + sk * cp ) * v[0] +
              ( -sk * ct * sp + ck * cp ) * v[1] + st * sp * v[2];
    temp[2] = -ck * st * v[0] + sk * st * v[1] + ct * v[2];

    v[0] = temp[0];
    v[1] = temp[1];
    v[2] = temp[2];
}

ostream& operator<<( ostream& s, const EvtVector3C& v )
{
    s << "(" << v.v[0] << "," << v.v[1] << "," << v.v[2] << ")";

    return s;
}
