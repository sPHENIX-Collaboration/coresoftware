
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

#include "EvtGenExternal/EvtPythia.hh"

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinDensity.hh"

#include "EvtGenModels/EvtAbsExternalGen.hh"

#include "EvtGenExternal/EvtExternalGenFactory.hh"

#include <cmath>
#include <iostream>

EvtPythia::EvtPythia()
{
    // Set the Pythia engine to a null pointer at first.
    // When we do the decay, we retrieve the pointer to the Pythia engine
    // and use that for all decays. All clones will use the same Pythia engine.
    _pythiaEngine = 0;
}

EvtPythia::~EvtPythia()
{
    _commandList.clear();
}

std::string EvtPythia::getName()
{
    return "PYTHIA";
}

EvtDecayBase* EvtPythia::clone()
{
    return new EvtPythia();
}

void EvtPythia::init()
{
    // Do not check for any arguments. The PythiaEngine will check
    // to see if there is an integer specifying the decay physics,
    // otherwise it just uses phase-space.
}

void EvtPythia::initProbMax()
{
    noProbMax();
}

void EvtPythia::decay( EvtParticle* p )
{
    // We have to initialise the Pythia engine after the decay.dec files have been read in,
    // since we will be modifying Pythia data tables, and that is only possible once we have
    // defined all Pythia-type decays we want to use.
    // We check to see if the engine has been created before doing the decay.
    // This should only create the full Pythia engine once, and all clones will point to the same engine.

    if ( !_pythiaEngine ) {
        _pythiaEngine = EvtExternalGenFactory::getInstance()->getGenerator(
            EvtExternalGenFactory::PythiaGenId );
    }

    if ( _pythiaEngine ) {
        _pythiaEngine->doDecay( p );
    }

    this->fixPolarisations( p );
}

void EvtPythia::fixPolarisations( EvtParticle* p )
{
    // Special case to handle the J/psi polarisation

    if ( !p ) {
        return;
    }

    int nDaug = p->getNDaug();
    int i( 0 );

    static EvtId Jpsi = EvtPDL::getId( "J/psi" );

    for ( i = 0; i < nDaug; i++ ) {
        EvtParticle* theDaug = p->getDaug( i );

        if ( theDaug ) {
            if ( theDaug->getId() == Jpsi ) {
                EvtSpinDensity rho;

                rho.setDim( 3 );
                rho.set( 0, 0, 0.5 );
                rho.set( 0, 1, 0.0 );
                rho.set( 0, 2, 0.0 );

                rho.set( 1, 0, 0.0 );
                rho.set( 1, 1, 1.0 );
                rho.set( 1, 2, 0.0 );

                rho.set( 2, 0, 0.0 );
                rho.set( 2, 1, 0.0 );
                rho.set( 2, 2, 0.5 );

                EvtVector4R p4Psi = theDaug->getP4();

                double alpha = atan2( p4Psi.get( 2 ), p4Psi.get( 1 ) );
                double beta = acos( p4Psi.get( 3 ) / p4Psi.d3mag() );

                theDaug->setSpinDensityForwardHelicityBasis( rho, alpha, beta,
                                                             0.0 );
                setDaughterSpinDensity( i );
            }
        }
    }
}

std::string EvtPythia::commandName()
{
    // Allow backward compatibility for decay.dec files
    // having JetSetPar parameters. They are obsolete for Pythia 8,
    // since the JetSet-type array variables do not exist.
    // Need to think about including user defined parameters in
    // EvtPythiaEngine::updatePhysicsParameters().
    return std::string( "JetSetPar" );
}

void EvtPythia::command( std::string cmd )
{
    // Locally store commands in a vector
    _commandList.push_back( cmd );
}
