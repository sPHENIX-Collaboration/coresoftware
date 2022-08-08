
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

#include "EvtGenBase/EvtScalarParticle.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iostream>
#include <math.h>

void EvtScalarParticle::init( EvtId part_n, double e, double px, double py,
                              double pz )
{
    _validP4 = true;
    setp( e, px, py, pz );
    setpart_num( part_n );

    setLifetime();
}

void EvtScalarParticle::init( EvtId part_n, const EvtVector4R& p4 )
{
    _validP4 = true;
    setp( p4 );
    setpart_num( part_n );

    setLifetime();
}

EvtSpinDensity EvtScalarParticle::rotateToHelicityBasis() const
{
    EvtSpinDensity R;
    R.setDim( 1 );

    R.set( 0, 0, 1.0 );

    return R;
}

EvtSpinDensity EvtScalarParticle::rotateToHelicityBasis( double, double,
                                                         double ) const
{
    EvtSpinDensity R;
    R.setDim( 1 );

    R.set( 0, 0, 1.0 );

    return R;
}
