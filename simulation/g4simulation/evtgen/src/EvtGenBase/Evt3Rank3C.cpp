
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

#include "EvtGenBase/Evt3Rank3C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor3C.hh"
#include "EvtGenBase/EvtVector3C.hh"

#include <iostream>
#include <math.h>

Evt3Rank3C::Evt3Rank3C( const Evt3Rank3C& t1 )
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] = t1.t[i][j][k];
            }
        }
    }
}

Evt3Rank3C& Evt3Rank3C::operator=( const Evt3Rank3C& t1 )
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] = t1.t[i][j][k];
            }
        }
    }
    return *this;
}

Evt3Rank3C Evt3Rank3C::conj() const
{
    Evt3Rank3C temp;

    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                temp.set( i, j, k, ::conj( t[i][j][k] ) );
            }
        }
    }
    return temp;
}

void Evt3Rank3C::zero()
{
    int i, j, k;
    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] = EvtComplex( 0.0, 0.0 );
            }
        }
    }
}

Evt3Rank3C::Evt3Rank3C()
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] = EvtComplex( 0.0, 0.0 );
            }
        }
    }
}

std::ostream& operator<<( std::ostream& s, const Evt3Rank3C& t2 )
{
    int i, j, k;
    for ( k = 0; k < 3; k++ ) {
        for ( i = 0; i < 3; i++ ) {
            for ( j = 0; j < 3; j++ ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << t2.t[k][i][j];
            }
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "\n";
        }
    }
    return s;
}

Evt3Rank3C& Evt3Rank3C::operator+=( const Evt3Rank3C& t2 )
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] += t2.t[i][j][k];
            }
        }
    }
    return *this;
}

Evt3Rank3C& Evt3Rank3C::operator-=( const Evt3Rank3C& t2 )
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] -= t2.t[i][j][k];
            }
        }
    }

    return *this;
}

Evt3Rank3C& Evt3Rank3C::operator*=( const EvtComplex& c )
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] *= c;
            }
        }
    }
    return *this;
}

Evt3Rank3C& Evt3Rank3C::operator*=( const double c )
{
    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                t[i][j][k] *= c;
            }
        }
    }

    return *this;
}

Evt3Rank3C conj( const Evt3Rank3C& t2 )
{
    Evt3Rank3C temp;

    int i, j, k;

    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 3; k++ ) {
                temp.t[i][j][k] = ::conj( t2.t[i][j][k] );
            }
        }
    }
    return temp;
}

EvtTensor3C Evt3Rank3C::cont1( const EvtVector3C& v ) const
{
    EvtTensor3C temp;

    int i, k;

    for ( i = 0; i < 3; i++ ) {
        for ( k = 0; k < 3; k++ ) {
            temp.set( i, k,
                      t[0][i][k] * v.get( 0 ) + t[1][i][k] * v.get( 1 ) +
                          t[2][i][k] * v.get( 2 ) );
        }
    }
    return temp;
}

EvtTensor3C Evt3Rank3C::cont2( const EvtVector3C& v ) const
{
    EvtTensor3C temp;

    int i, k;

    for ( i = 0; i < 3; i++ ) {
        for ( k = 0; k < 3; k++ ) {
            temp.set( i, k,
                      t[i][0][k] * v.get( 0 ) + t[i][1][k] * v.get( 1 ) +
                          t[i][2][k] * v.get( 2 ) );
        }
    }
    return temp;
}

EvtTensor3C Evt3Rank3C::cont3( const EvtVector3C& v ) const
{
    EvtTensor3C temp;

    int i, k;

    for ( i = 0; i < 3; i++ ) {
        for ( k = 0; k < 3; k++ ) {
            temp.set( i, k,
                      t[i][k][0] * v.get( 0 ) + t[i][k][1] * v.get( 1 ) +
                          t[i][k][2] * v.get( 2 ) );
        }
    }
    return temp;
}

EvtTensor3C Evt3Rank3C::cont1( const EvtVector3R& v ) const
{
    EvtTensor3C temp;

    int i, k;

    for ( i = 0; i < 3; i++ ) {
        for ( k = 0; k < 3; k++ ) {
            temp.set( i, k,
                      t[0][i][k] * v.get( 0 ) + t[1][i][k] * v.get( 1 ) +
                          t[2][i][k] * v.get( 2 ) );
        }
    }
    return temp;
}

EvtTensor3C Evt3Rank3C::cont2( const EvtVector3R& v ) const
{
    EvtTensor3C temp;

    int i, k;

    for ( i = 0; i < 3; i++ ) {
        for ( k = 0; k < 3; k++ ) {
            temp.set( i, k,
                      t[i][0][k] * v.get( 0 ) + t[i][1][k] * v.get( 1 ) +
                          t[i][2][k] * v.get( 2 ) );
        }
    }
    return temp;
}

EvtTensor3C Evt3Rank3C::cont3( const EvtVector3R& v ) const
{
    EvtTensor3C temp;

    int i, k;

    for ( i = 0; i < 3; i++ ) {
        for ( k = 0; k < 3; k++ ) {
            temp.set( i, k,
                      t[i][k][0] * v.get( 0 ) + t[i][k][1] * v.get( 1 ) +
                          t[i][k][2] * v.get( 2 ) );
        }
    }
    return temp;
}
