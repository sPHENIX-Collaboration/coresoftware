
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

#include "EvtGenModels/EvtItgAbsIntegrator.hh"

#include "EvtGenBase/EvtPatches.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtItgAbsFunction.hh"

#include <iostream>
#include <math.h>
using std::endl;

EvtItgAbsIntegrator::EvtItgAbsIntegrator( const EvtItgAbsFunction& theFunction ) :
    _myFunction( theFunction )
{
}

double EvtItgAbsIntegrator::normalisation() const
{
    return evaluateIt( _myFunction.lowerRange(), _myFunction.upperRange() );
}

double EvtItgAbsIntegrator::evaluate( double lower, double upper ) const
{
    double newLower( lower ), newUpper( upper );

    boundsCheck( newLower, newUpper );

    return evaluateIt( newLower, newUpper );
}

double EvtItgAbsIntegrator::trapezoid( double lower, double higher, int n,
                                       double& result ) const
{
    if ( n == 1 )
        return 0.5 * ( higher - lower ) *
               ( _myFunction( lower ) + _myFunction( higher ) );

    int it, j;

    for ( it = 1, j = 1; j < n - 1; j++ )
        it <<= 1;

    double itDouble( it );

    double sum( 0.0 );

    double deltaX( ( higher - lower ) / itDouble );

    double x( lower + 0.5 * deltaX );

    for ( j = 1; j <= it; j++ ) {
        sum += _myFunction( x );
        x += deltaX;
    }

    result = 0.5 * ( result + ( higher - lower ) * sum / itDouble );

    return result;
}

void EvtItgAbsIntegrator::boundsCheck( double& lower, double& upper ) const
{
    if ( lower < _myFunction.lowerRange() ) {
        EvtGenReport( EVTGEN_WARNING, "EvtGen" )
            << "Warning in EvtItgAbsIntegrator::evaluate.  Lower bound "
            << lower << " of integral "
            << " is less than lower bound " << _myFunction.lowerRange()
            << " of function.  No contribution from this range will be counted."
            << endl;
        lower = _myFunction.lowerRange();
    }

    if ( upper > _myFunction.upperRange() ) {
        EvtGenReport( EVTGEN_WARNING, "EvtGen" )
            << "Warning in EvtItgAbsIntegrator::evaluate.  Upper bound "
            << upper << " of integral "
            << " is greater than upper bound " << _myFunction.upperRange()
            << " of function.  No contribution from this range will be counted."
            << endl;
        upper = _myFunction.upperRange();
    }
}
