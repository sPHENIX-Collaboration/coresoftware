
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

#include "EvtGenModels/EvtPhiDalitz.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <math.h>
#include <stdlib.h>
#include <string>

// Implementation of KLOE measurement
// PL B561: 55-60 (2003) + Erratum B609:449-450 (2005)
// or hep-ex/0303016v2

std::string EvtPhiDalitz::getName()
{
    return "PHI_DALITZ";
}

EvtDecayBase* EvtPhiDalitz::clone()
{
    return new EvtPhiDalitz;
}

void EvtPhiDalitz::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::VECTOR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );

    _mRho = 0.7758;
    _gRho = 0.1439;
    _aD = 0.78;
    _phiD = -2.47;
    _aOmega = 0.0071;
    _phiOmega = -0.22;

    _locPip = -1;
    _locPim = -1;
    _locPi0 = -1;

    for ( int i = 0; i < 3; i++ ) {
        if ( getDaug( i ) == EvtPDL::getId( "pi+" ) )
            _locPip = i;
        if ( getDaug( i ) == EvtPDL::getId( "pi-" ) )
            _locPim = i;
        if ( getDaug( i ) == EvtPDL::getId( "pi0" ) )
            _locPi0 = i;
    }
    if ( _locPip == -1 || _locPim == -1 || _locPi0 == -1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << getModelName()
            << "generator expects daughters to be pi+ pi- pi0\n";
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Found " << EvtPDL::name( getDaug( 0 ) ) << " "
            << EvtPDL::name( getDaug( 1 ) ) << " "
            << EvtPDL::name( getDaug( 2 ) ) << std::endl;
    }
}

void EvtPhiDalitz::decay( EvtParticle* p )
{
    EvtId PIP = EvtPDL::getId( "pi+" );
    EvtId PIM = EvtPDL::getId( "pi-" );
    EvtId PIZ = EvtPDL::getId( "pi0" );
    EvtId OMEGA = EvtPDL::getId( "omega" );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R Ppip = p->getDaug( _locPip )->getP4();
    EvtVector4R Ppim = p->getDaug( _locPim )->getP4();
    EvtVector4R Ppi0 = p->getDaug( _locPi0 )->getP4();
    EvtVector4R Qp = ( Ppim + Ppi0 );
    EvtVector4R Qm = ( Ppip + Ppi0 );
    EvtVector4R Q0 = ( Ppip + Ppim );
    double m2_pip = pow( EvtPDL::getMeanMass( PIP ), 2 );
    double m2_pim = pow( EvtPDL::getMeanMass( PIM ), 2 );
    double m2_pi0 = pow( EvtPDL::getMeanMass( PIZ ), 2 );
    double M2rhop = pow( _mRho, 2 );
    double M2rhom = pow( _mRho, 2 );
    double M2rho0 = pow( _mRho, 2 );
    double M2omega = pow( EvtPDL::getMeanMass( OMEGA ), 2 );

    double Wrhop = _gRho;
    double Wrhom = _gRho;
    double Wrho0 = _gRho;
    double Womega = EvtPDL::getWidth( OMEGA );

    EvtComplex Atot( 0, 0 );

    //Rho+ Risonance Amplitude
    double Gp = Wrhop *
                pow( ( ( Qp.mass2() - m2_pim - m2_pi0 ) / 2 - M2rhop / 4 ) /
                         ( M2rhop / 4 - ( m2_pim + m2_pi0 ) / 2 ),
                     3 / 2 ) *
                ( M2rhop / Qp.mass2() );
    EvtComplex Drhop( ( Qp.mass2() - M2rhop ), Qp.mass() * Gp );
    EvtComplex A1( M2rhop / Drhop );

    //Rho- Risonance Amplitude
    double Gm = Wrhom *
                pow( ( ( Qm.mass2() - m2_pip - m2_pi0 ) / 2 - M2rhom / 4 ) /
                         ( M2rhom / 4 - ( m2_pip + m2_pi0 ) / 2 ),
                     3 / 2 ) *
                ( M2rhom / Qm.mass2() );
    EvtComplex Drhom( ( Qm.mass2() - M2rhom ), Qm.mass() * Gm );
    EvtComplex A2( M2rhom / Drhom );

    //Rho0 Risonance Amplitude
    double G0 = Wrho0 *
                pow( ( ( Q0.mass2() - m2_pip - m2_pim ) / 2 - M2rho0 / 4 ) /
                         ( M2rho0 / 4 - ( m2_pip + m2_pim ) / 2 ),
                     3 / 2 ) *
                ( M2rho0 / Q0.mass2() );
    EvtComplex Drho0( ( Q0.mass2() - M2rho0 ), Q0.mass() * G0 );
    EvtComplex A3( M2rho0 / Drho0 );

    //Omega Risonance Amplitude
    EvtComplex OmegaPhase( 0, _phiOmega );
    EvtComplex DOmega( ( Q0.mass2() - M2omega ), Q0.mass() * Womega );
    EvtComplex A4( _aOmega * M2omega * exp( OmegaPhase ) / DOmega );

    //Direct Decay Amplitude
    EvtComplex DirPhase( 0, _phiD );
    EvtComplex A5( _aD * exp( DirPhase ) );

    Atot = A1 + A2 + A3 + A4 + A5;

    vertex( 0, Atot );
    vertex( 1, Atot );
    vertex( 2, Atot );

    return;
}
