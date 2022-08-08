
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

#include "EvtGenBase/EvtPto3PAmpSmpResolution.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPto3PAmp.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using EvtCyclic3::Index;
using EvtCyclic3::Pair;
using std::cout;
using std::endl;

EvtPto3PAmpSmpResolution::EvtPto3PAmpSmpResolution( EvtDalitzPlot dp,
                                                    Pair pairAng, Pair pairRes,
                                                    EvtSpinType::spintype spin,
                                                    const EvtPropagator& prop,
                                                    NumType typeN ) :
    EvtPto3PAmp( dp, pairAng, pairRes, spin, prop, typeN )
{
}

EvtComplex EvtPto3PAmpSmpResolution::evalPropagator( double m ) const
{
    EvtComplex prop( 0, 0 );

    if ( _sigma > 0 ) {    // convolved
        int nconv = 20;
        double min = m + _bias - _sigma * 2.5;
        double max = m + _bias + _sigma * 2.5;
        double dm = ( max - min ) / nconv;
        static double sqrt2pi = sqrt( 2 * 3.14159 );
        double ifact = 1. / ( sqrt2pi * _sigma );
        for ( int i = 0; i < nconv; i++ ) {
            double mprime = min + dm * ( i + 0.5 );
            double t = ( mprime - m ) / _sigma;
            prop += ifact * exp( -0.5 * t * t ) *
                    EvtPto3PAmp::evalPropagator( m ) * dm;
        }
    } else {
        prop = EvtPto3PAmp::evalPropagator( m );
    }

    return prop;
}
