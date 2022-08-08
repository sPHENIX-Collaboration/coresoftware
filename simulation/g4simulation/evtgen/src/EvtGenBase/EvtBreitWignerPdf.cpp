
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

#include "EvtGenBase/EvtBreitWignerPdf.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <math.h>
#include <stdio.h>

EvtBreitWignerPdf::EvtBreitWignerPdf( double min, double max, double m0,
                                      double g0 ) :
    EvtIntegPdf1D( min, max ), _m0( m0 ), _g0( g0 )
{
}

EvtBreitWignerPdf::EvtBreitWignerPdf( const EvtBreitWignerPdf& other ) :
    EvtIntegPdf1D( other ), _m0( other._m0 ), _g0( other._g0 )
{
}

double EvtBreitWignerPdf::pdf( const EvtPoint1D& x ) const
{
    double m = x.value();
    if ( ( 0 == ( m - _m0 ) ) && ( 0. == _g0 ) ) {
        printf( "Delta function Breit-Wigner\n" );
        assert( 0 );
    }

    double ret = _g0 / EvtConst::twoPi /
                 ( ( m - _m0 ) * ( m - _m0 ) + _g0 * _g0 / 4 );

    return ret;
}

double EvtBreitWignerPdf::pdfIntegral( double m ) const
{
    double itg = 0;
    if ( _g0 == 0 ) {
        if ( m > _m0 )
            itg = 1.;
        else if ( m < _m0 )
            itg = 0.;
        else
            itg = 0.5;
    } else
        itg = atan( ( m - _m0 ) / ( _g0 / 2. ) ) / EvtConst::pi + 0.5;

    return itg;
}

double EvtBreitWignerPdf::pdfIntegralInverse( double x ) const
{
    if ( x < 0 || x > 1 ) {
        printf( "Invalid integral value %f\n", x );
        assert( 0 );
    }

    double m = _m0;
    if ( _g0 != 0 )
        m = _m0 + ( _g0 / 2. ) * tan( EvtConst::pi * ( x - 0.5 ) );

    return m;
}
