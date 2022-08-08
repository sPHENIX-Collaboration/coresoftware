
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

#include "EvtGenBase/EvtHighSpinParticle.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <assert.h>
#include <iostream>
#include <math.h>

void EvtHighSpinParticle::init( EvtId id, const EvtVector4R& p4 )
{
    _validP4 = true;
    setp( p4 );
    setpart_num( id );

    setLifetime();
}

EvtSpinDensity EvtHighSpinParticle::rotateToHelicityBasis() const
{
    int n = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getId() ) );

    EvtSpinDensity R;
    R.setDiag( n );

    return R;
}

EvtSpinDensity EvtHighSpinParticle::rotateToHelicityBasis( double alpha,
                                                           double beta,
                                                           double gamma ) const
{
    int i, j;

    int n = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getId() ) );

    EvtSpinDensity R;

    R.setDim( n );

    int J2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getId() ) );

    assert( n == J2 + 1 );

    auto lambda2 = [J2]( int i ) { return J2 - i * 2; };

    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < n; j++ ) {
            R.set( i, j,
                   EvtdFunction::d( J2, lambda2( j ), lambda2( i ), beta ) *
                       exp( EvtComplex( 0.0, 0.5 * ( alpha * lambda2( i ) -
                                                     gamma * lambda2( j ) ) ) ) );
        }
    }

    return R;
}
