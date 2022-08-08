
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

#include "EvtGenBase/EvtPropBreitWignerRel.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <math.h>

EvtPropBreitWignerRel::EvtPropBreitWignerRel( double m0, double g0 ) :
    EvtPropagator( m0, g0 )
{
}

EvtAmplitude<EvtPoint1D>* EvtPropBreitWignerRel::clone() const
{
    return new EvtPropBreitWignerRel( *this );
}

EvtComplex EvtPropBreitWignerRel::amplitude( const EvtPoint1D& x ) const
{
    double m = x.value();
    return 1. / ( _m0 * _m0 - m * m - EvtComplex( 0., _m0 * _g0 ) );
}
