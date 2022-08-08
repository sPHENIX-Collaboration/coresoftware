
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

#include "EvtGenBase/EvtTensor3C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector3C.hh"

#include <iostream>
#include <math.h>
using std::endl;
using std::ostream;

EvtTensor3C::EvtTensor3C( const EvtTensor3C& t1 )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] = t1.t[i][j];
        }
    }
}

EvtTensor3C::EvtTensor3C( double d11, double d22, double d33 )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] = 0.0;
        }
    }

    t[0][0] = d11;
    t[1][1] = d22;
    t[2][2] = d33;
}

EvtTensor3C& EvtTensor3C::operator=( const EvtTensor3C& t1 )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] = t1.t[i][j];
        }
    }
    return *this;
}

EvtTensor3C EvtTensor3C::conj() const
{
    EvtTensor3C temp;

    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            temp.set( j, i, ::conj( t[i][j] ) );
        }
    }
    return temp;
}

void EvtTensor3C::zero()
{
    int i, j;
    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] = EvtComplex( 0.0, 0.0 );
        }
    }
}

EvtTensor3C::EvtTensor3C()
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] = EvtComplex( 0.0, 0.0 );
        }
    }
}

EvtTensor3C EvtTensor3C::operator+=( const EvtTensor3C& t2 )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] += t2.t[i][j];
        }
    }
    return *this;
}

EvtTensor3C EvtTensor3C::operator-=( const EvtTensor3C& t2 )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] -= t2.t[i][j];
        }
    }
    return *this;
}

EvtTensor3C EvtTensor3C::operator*=( const EvtComplex& c )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] *= c;
        }
    }
    return *this;
}

EvtTensor3C EvtTensor3C::operator*=( const double c )
{
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] *= EvtComplex( c );
        }
    }
    return *this;
}

EvtTensor3C EvtGenFunctions::directProd( const EvtVector3C& c1,
                                         const EvtVector3C& c2 )
{
    EvtTensor3C temp;
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            temp.set( i, j, c1.get( i ) * c2.get( j ) );
        }
    }
    return temp;
}

EvtTensor3C EvtGenFunctions::directProd( const EvtVector3C& c1,
                                         const EvtVector3R& c2 )
{
    EvtTensor3C temp;
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            temp.set( i, j, c1.get( i ) * c2.get( j ) );
        }
    }
    return temp;
}

EvtTensor3C EvtGenFunctions::directProd( const EvtVector3R& c1,
                                         const EvtVector3R& c2 )
{
    EvtTensor3C temp;
    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            temp.t[i][j] = EvtComplex( c1.get( i ) * c2.get( j ), 0.0 );
        }
    }
    return temp;
}

EvtTensor3C conj( const EvtTensor3C& t2 )
{
    EvtTensor3C temp;

    int i, j;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            temp.set( i, j, ::conj( ( t2.get( i, j ) ) ) );
        }
    }

    return temp;
}

EvtTensor3C cont22( const EvtTensor3C& t1, const EvtTensor3C& t2 )
{
    EvtTensor3C temp;

    int i, j;
    EvtComplex c;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            c = t1.get( i, 0 ) * t2.get( j, 0 ) +
                t1.get( i, 1 ) * t2.get( j, 1 ) + t1.get( i, 2 ) * t2.get( j, 2 );
            temp.set( i, j, c );
        }
    }

    return temp;
}

EvtTensor3C cont11( const EvtTensor3C& t1, const EvtTensor3C& t2 )
{
    EvtTensor3C temp;

    int i, j;
    EvtComplex c;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            c = t1.get( 0, i ) * t2.get( 0, j ) +
                t1.get( 1, i ) * t2.get( 1, j ) + t1.get( 2, i ) * t2.get( 2, j );
            temp.set( i, j, c );
        }
    }

    return temp;
}

EvtVector3C EvtTensor3C::cont1( const EvtVector3C& v ) const
{
    EvtVector3C temp;

    int i;

    for ( i = 0; i < 3; i++ ) {
        temp.set( i, t[0][i] * v.get( 0 ) + t[1][i] * v.get( 1 ) +
                         t[2][i] * v.get( 2 ) );
    }

    return temp;
}

EvtVector3C EvtTensor3C::cont2( const EvtVector3C& v ) const
{
    EvtVector3C temp;

    int i;

    for ( i = 0; i < 3; i++ ) {
        temp.set( i, t[i][0] * v.get( 0 ) + t[i][1] * v.get( 1 ) +
                         t[i][2] * v.get( 2 ) );
    }

    return temp;
}

EvtVector3C EvtTensor3C::cont1( const EvtVector3R& v ) const
{
    EvtVector3C temp;

    int i;

    for ( i = 0; i < 3; i++ ) {
        temp.set( i, t[0][i] * v.get( 0 ) + t[1][i] * v.get( 1 ) +
                         t[2][i] * v.get( 2 ) );
    }

    return temp;
}

EvtVector3C EvtTensor3C::cont2( const EvtVector3R& v ) const
{
    EvtVector3C temp;

    int i;

    for ( i = 0; i < 3; i++ ) {
        temp.set( i, t[i][0] * v.get( 0 ) + t[i][1] * v.get( 1 ) +
                         t[i][2] * v.get( 2 ) );
    }

    return temp;
}

EvtTensor3C EvtGenFunctions::eps( const EvtVector3R& v )
{
    EvtTensor3C temp;

    temp.t[0][0] = 0.0;
    temp.t[1][1] = 0.0;
    temp.t[2][2] = 0.0;

    temp.t[0][1] = v.get( 2 );
    temp.t[0][2] = -v.get( 1 );

    temp.t[1][0] = -v.get( 2 );
    temp.t[1][2] = v.get( 0 );

    temp.t[2][0] = v.get( 1 );
    temp.t[2][1] = -v.get( 0 );

    return temp;
}

const EvtTensor3C& EvtTensor3C::id()
{
    static EvtTensor3C identity( 1.0, 1.0, 1.0 );

    return identity;
}

ostream& operator<<( ostream& s, const EvtTensor3C& v )
{
    s << endl
      << "(" << v.t[0][0] << "," << v.t[0][1] << "," << v.t[0][2] << ")";
    s << endl
      << "(" << v.t[1][0] << "," << v.t[1][1] << "," << v.t[1][2] << ")";
    s << endl
      << "(" << v.t[2][0] << "," << v.t[2][1] << "," << v.t[2][2] << ")" << endl;

    return s;
}

EvtTensor3C EvtGenFunctions::rotateEuler( const EvtTensor3C& v, double alpha,
                                          double beta, double gamma )
{
    EvtTensor3C tmp( v );
    tmp.applyRotateEuler( alpha, beta, gamma );
    return tmp;
}

void EvtTensor3C::applyRotateEuler( double phi, double theta, double ksi )
{
    EvtComplex temp[3][3];
    double sp, st, sk, cp, ct, ck;
    double r[3][3];
    int i, j, k;

    sp = sin( phi );
    st = sin( theta );
    sk = sin( ksi );
    cp = cos( phi );
    ct = cos( theta );
    ck = cos( ksi );

    r[0][0] = ck * ct * cp - sk * sp;
    r[0][1] = ck * ct * sp + sk * cp;
    r[0][2] = -ck * st;

    r[1][0] = -sk * ct * cp - ck * sp;
    r[1][1] = -sk * ct * sp + ck * cp;
    r[1][2] = sk * st;

    r[2][0] = st * cp;
    r[2][1] = st * sp;
    r[2][2] = ct;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            temp[i][j] = 0.0;
            for ( k = 0; k < 3; k++ ) {
                temp[i][j] += r[i][k] * t[k][j];
            }
        }
    }

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            t[i][j] = 0.0;
            for ( k = 0; k < 3; k++ ) {
                t[i][j] += r[i][k] * temp[j][k];
            }
        }
    }
}
