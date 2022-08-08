
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

#include "EvtGenBase/EvtDalitzResPdf.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"

#include <math.h>
#include <stdio.h>
using namespace EvtCyclic3;

EvtDalitzResPdf::EvtDalitzResPdf( const EvtDalitzPlot& dp, double _m0,
                                  double _g0, EvtCyclic3::Pair pair ) :
    EvtPdf<EvtDalitzPoint>(), _dp( dp ), _m0( _m0 ), _g0( _g0 ), _pair( pair )
{
}

EvtValError EvtDalitzResPdf::compute_integral( int N ) const
{
    assert( N != 0 );

    EvtCyclic3::Pair i = _pair;
    EvtCyclic3::Pair j = EvtCyclic3::next( i );

    // Trapezoidal integral

    double dh = ( _dp.qAbsMax( j ) - _dp.qAbsMin( j ) ) / ( (double)N );
    double sum = 0;

    int ii;
    for ( ii = 1; ii < N; ii++ ) {
        double x = _dp.qAbsMin( j ) + ii * dh;
        double min = ( _dp.qMin( i, j, x ) - _m0 * _m0 ) / _m0 / _g0;
        double max = ( _dp.qMax( i, j, x ) - _m0 * _m0 ) / _m0 / _g0;
        double itg = 1 / EvtConst::pi * ( atan( max ) - atan( min ) );
        sum += itg;
    }
    EvtValError ret( sum * dh, 0. );

    return ret;
}

EvtDalitzPoint EvtDalitzResPdf::randomPoint()
{
    // Random point generation must be done in a box encompassing the
    // Dalitz plot

    EvtCyclic3::Pair i = _pair;
    EvtCyclic3::Pair j = EvtCyclic3::next( i );
    double min = 1 / EvtConst::pi *
                 atan( ( _dp.qAbsMin( i ) - _m0 * _m0 ) / _m0 / _g0 );
    double max = 1 / EvtConst::pi *
                 atan( ( _dp.qAbsMax( i ) - _m0 * _m0 ) / _m0 / _g0 );

    int n = 0;
    while ( n++ < 1000 ) {
        double qj = EvtRandom::Flat( _dp.qAbsMin( j ), _dp.qAbsMax( j ) );
        double r = EvtRandom::Flat( min, max );
        double qi = tan( EvtConst::pi * r ) * _g0 * _m0 + _m0 * _m0;
        EvtDalitzCoord x( i, qi, j, qj );
        EvtDalitzPoint ret( _dp, x );
        if ( ret.isValid() )
            return ret;
    }

    // All generated points turned out to be outside of the Dalitz plot
    // (in the outer box)

    printf( "No point generated for dalitz plot after 1000 tries\n" );
    return EvtDalitzPoint( 0., 0., 0., 0., 0., 0. );
}

double EvtDalitzResPdf::pdf( const EvtDalitzPoint& x ) const
{
    EvtCyclic3::Pair i = _pair;
    double dq = x.q( i ) - _m0 * _m0;
    return 1 / EvtConst::pi * _g0 * _m0 / ( dq * dq + _g0 * _g0 * _m0 * _m0 );
}

double EvtDalitzResPdf::pdfMaxValue() const
{
    return 1 / ( EvtConst::pi * _g0 * _m0 );
}
