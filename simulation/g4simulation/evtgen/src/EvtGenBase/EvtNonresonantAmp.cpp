
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

#include "EvtGenBase/EvtNonresonantAmp.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using EvtCyclic3::Index;
using EvtCyclic3::Pair;
using std::endl;

EvtNonresonantAmp::EvtNonresonantAmp( EvtDalitzPlot* dp,
                                      EvtPto3PAmp::NumType type,
                                      EvtCyclic3::Pair pair1, double par1,
                                      EvtCyclic3::Pair pair2, double par2,
                                      EvtSpinType::spintype spin ) :
    EvtAmplitude<EvtDalitzPoint>(),
    _dalitzSpace{dp},
    _type( type ),
    _pair1( pair1 ),
    _pair2( pair2 ),
    _par1( par1 ),
    _par2( par2 ),
    _spin( spin )
{
}

EvtComplex EvtNonresonantAmp::amplitude( const EvtDalitzPoint& dalitzPoint ) const
{
    // flat model
    if ( _type == EvtPto3PAmp::NONRES ) {
        return 1;
    }

    // "linear model" (prop. to m^2)
    else if ( _type == EvtPto3PAmp::NONRES_LIN ) {
        return dalitzPoint.q( _pair1 );
    }

    // Chen-Chua-Soni
    else if ( _type == EvtPto3PAmp::NONRES_CCS ) {
        double s = dalitzPoint.q( _pair1 );
        double smin = _dalitzSpace->qAbsMin( _pair1 );
        return sqrt( s - smin ) / ( s * log( s * _par1 ) );
    }

    // exp{par*m^2) (Belle model, Garmash et al, PRD71)
    else if ( _type == EvtPto3PAmp::NONRES_EXP ) {
        return exp( _par1 * dalitzPoint.q( _pair1 ) );
    }

    // exp(par1*m12^2 + par2*m13^2) (Belle model, Garmash et al, PRD71)
    else if ( _type == EvtPto3PAmp::NONRES_EXP_ADD ) {
        return exp( _par1 * dalitzPoint.q( _pair1 ) +
                    _par2 * dalitzPoint.q( _pair2 ) );
    }

    // Laura model (P.Harrison et al, BAD806)
    else if ( _type == EvtPto3PAmp::NONRES_LAURA ) {
        double m = sqrt( dalitzPoint.q( _pair1 ) );
        double mmin = sqrt( _dalitzSpace->qAbsMin( _pair1 ) );
        double dm = m - mmin;
        assert( dm > 0 );
        double cosTh = 1;
        int ispin = EvtSpinType::getSpin2( _spin );
        if ( ispin > 0 ) {
            cosTh = dalitzPoint.cosTh( EvtCyclic3::next( _pair1 ), _pair1 );
            if ( ispin > 2 )
                cosTh *= cosTh;
        }
        return pow( dm, _par1 ) * exp( dm * _par2 ) * cosTh;
    }

    return 0;
}
