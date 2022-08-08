
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

#include "EvtGenBase/EvtDiLog.hh"

#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : EvtDiLog
//
// 2007-01-23 : Patrick Robbe
//-----------------------------------------------------------------------------

double EvtDiLog::DiLog( double x )
{
    double h, t, y, s, a, alfa, b0, b1, b2;
    if ( x == 1. )
        h = PI6;
    else if ( x == -1. )
        h = -PI12;
    else {
        t = -x;
        if ( t <= -2. ) {
            y = -1. / ( 1. + t );
            s = 1.;
            a = -PI3 + HF * ( std::pow( log( -t ), 2 ) -
                              std::pow( log( 1. + 1. / t ), 2 ) );
        } else if ( t < -1. ) {
            y = -1. - t;
            s = -1.;
            a = log( -t );
            a = -PI6 + a * ( a + log( 1. + 1. / t ) );
        } else if ( t <= -HF ) {
            y = -( 1. + t ) / t;
            s = 1.;
            a = log( -t );
            a = -PI6 + a * ( -HF * a + log( 1. + t ) );
        } else if ( t < 0 ) {
            y = -t / ( 1. + t );
            s = -1.;
            a = HF * std::pow( log( 1. + t ), 2 );
        } else if ( t <= 1. ) {
            y = t;
            s = 1.;
            a = 0.;
        } else {
            y = 1. / t;
            s = -1.;
            a = PI6 + HF * std::pow( log( t ), 2 );
        }

        h = y + y - 1.;
        alfa = h + h;
        b1 = 0.;
        b2 = 0.;
        for ( int i = 19; i >= 0; --i ) {
            b0 = C[i] + alfa * b1 - b2;
            b2 = b1;
            b1 = b0;
        }

        h = -( s * ( b0 - h * b2 ) + a );
    }

    return h;
}
