
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

#include "EvtGenBase/EvtParticleFactory.hh"

#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtHighSpinParticle.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtNeutrinoParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtRaritaSchwingerParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtStringParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"

#include <sys/stat.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using std::endl;

EvtParticle* EvtParticleFactory::particleFactory( EvtSpinType::spintype spinType )
{
    if ( spinType == EvtSpinType::SCALAR ) {
        return new EvtScalarParticle;
    }

    if ( spinType == EvtSpinType::VECTOR ) {
        return new EvtVectorParticle;
    }
    if ( spinType == EvtSpinType::DIRAC ) {
        return new EvtDiracParticle;
    }
    if ( spinType == EvtSpinType::NEUTRINO ) {
        return new EvtNeutrinoParticle;
    }
    if ( spinType == EvtSpinType::PHOTON ) {
        return new EvtPhotonParticle;
    }
    if ( spinType == EvtSpinType::TENSOR ) {
        return new EvtTensorParticle;
    }
    if ( spinType == EvtSpinType::STRING ) {
        return new EvtStringParticle;
    }
    if ( spinType == EvtSpinType::RARITASCHWINGER ) {
        return new EvtRaritaSchwingerParticle;
    }
    if ( spinType == EvtSpinType::SPIN5HALF ) {
        return new EvtHighSpinParticle;
    }
    if ( spinType == EvtSpinType::SPIN3 ) {
        return new EvtHighSpinParticle;
    }
    if ( spinType == EvtSpinType::SPIN7HALF ) {
        return new EvtHighSpinParticle;
    }
    if ( spinType == EvtSpinType::SPIN4 ) {
        return new EvtHighSpinParticle;
    }

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Error in EvtParticleFactory::particleFactory" << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Tried to create non-existing particle"
        << " with spin type:" << spinType << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will terminate execution" << endl;

    ::abort();

    return 0;
}

EvtParticle* EvtParticleFactory::particleFactory( EvtId id, EvtVector4R p4,
                                                  EvtSpinDensity rho )
{
    EvtSpinType::spintype thisSpin = EvtPDL::getSpinType( id );

    if ( thisSpin == EvtSpinType::SCALAR ) {
        EvtScalarParticle* myPart;
        myPart = new EvtScalarParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }

    if ( thisSpin == EvtSpinType::VECTOR ) {
        EvtVectorParticle* myPart;
        myPart = new EvtVectorParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::DIRAC ) {
        EvtDiracParticle* myPart;
        myPart = new EvtDiracParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::NEUTRINO ) {
        EvtNeutrinoParticle* myPart;
        myPart = new EvtNeutrinoParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::PHOTON ) {
        EvtPhotonParticle* myPart;
        myPart = new EvtPhotonParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::TENSOR ) {
        EvtTensorParticle* myPart;
        myPart = new EvtTensorParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::STRING ) {
        EvtStringParticle* myPart;
        myPart = new EvtStringParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::SPIN3 ) {
        EvtHighSpinParticle* myPart;
        myPart = new EvtHighSpinParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::SPIN5HALF ) {
        EvtHighSpinParticle* myPart;
        myPart = new EvtHighSpinParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::SPIN7HALF ) {
        EvtHighSpinParticle* myPart;
        myPart = new EvtHighSpinParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::RARITASCHWINGER ) {
        EvtRaritaSchwingerParticle* myPart;
        myPart = new EvtRaritaSchwingerParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }
    if ( thisSpin == EvtSpinType::SPIN4 ) {
        EvtHighSpinParticle* myPart;
        myPart = new EvtHighSpinParticle;
        myPart->init( id, p4 );
        myPart->setSpinDensityForward( rho );
        return myPart;
    }

    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Error in EvtParticleFactory::particleFactory" << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Tried to create non-existing particle"
        << " with spin type:" << thisSpin
        << "  and name:" << EvtPDL::name( id ).c_str() << endl;
    EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will terminate execution" << endl;

    ::abort();

    return 0;
}

EvtParticle* EvtParticleFactory::particleFactory( EvtId id, EvtVector4R p4 )
{
    EvtSpinDensity rho;
    rho.setDiag( EvtSpinType::getSpinStates( EvtPDL::getSpinType( id ) ) );

    return particleFactory( id, p4, rho );
}
