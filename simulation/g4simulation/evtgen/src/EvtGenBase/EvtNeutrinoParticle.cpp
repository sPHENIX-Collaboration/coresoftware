
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

#include "EvtGenBase/EvtNeutrinoParticle.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iostream>
#include <math.h>
#include <stdlib.h>
using std::endl;

void EvtNeutrinoParticle::init( EvtId part_n, const EvtVector4R& p4 )
{
    _validP4 = true;
    setp( p4 );
    setpart_num( part_n );

    double e, px, py, pz;
    e = p4.get( 0 );
    px = p4.get( 1 );
    py = p4.get( 2 );
    pz = p4.get( 3 );

    if ( EvtPDL::getStdHep( part_n ) == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Error in EvtNeutrinoParticle::init, part_n=" << part_n.getId()
            << endl;
    }

    if ( EvtPDL::getStdHep( part_n ) > 0 ) {
        double beta, alpha, p2, norm;

        // See Sakurai p. 167-169
        // and Renton p. 126

        p2 = px * px + py * py + pz * pz;

        beta = acos( pz / sqrt( p2 ) );
        alpha = atan2( py, px );

        norm = sqrt( 2 * e );

        double cosb, sinb, cosa, sina;

        cosb = cos( 0.5 * beta );
        sinb = sin( 0.5 * beta );

        cosa = cos( 0.5 * alpha );
        sina = sin( 0.5 * alpha );

        spinor_parent.set( -norm * sinb * EvtComplex( cosa, -sina ),
                           norm * cosb * EvtComplex( cosa, sina ),
                           norm * sinb * EvtComplex( cosa, -sina ),
                           -norm * cosb * EvtComplex( cosa, sina ) );

    } else {
        px = -p4.get( 1 );
        py = -p4.get( 2 );
        pz = -p4.get( 3 );

        double pn, sqrpn;

        pn = e;
        sqrpn = sqrt( pn - pz );

        spinor_parent.set( ( 1.0 / sqrpn ) * EvtComplex( px, -py ),
                           EvtComplex( sqrpn, 0.0 ),
                           ( -1.0 / sqrpn ) * EvtComplex( px, -py ),
                           -EvtComplex( sqrpn, 0.0 ) );
    }

    setLifetime();
}

EvtDiracSpinor EvtNeutrinoParticle::spParentNeutrino() const
{
    return spinor_parent;
}

EvtDiracSpinor EvtNeutrinoParticle::spNeutrino() const
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Tried to get neutrino spinor in restframe";
    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will terminate execution.";

    ::abort();

    return spinor_rest;
}

EvtSpinDensity EvtNeutrinoParticle::rotateToHelicityBasis() const
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "rotateToHelicityBasis not implemented for neutrino.";
    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will terminate execution.";

    ::abort();

    EvtSpinDensity rho;
    return rho;
}

EvtSpinDensity EvtNeutrinoParticle::rotateToHelicityBasis( double, double,
                                                           double ) const
{
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "rotateToHelicityBasis(alpha,beta,gama) not implemented for neutrino.";
    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will terminate execution.";

    ::abort();

    EvtSpinDensity R;
    R.setDiag( 1 );

    return R;
}
