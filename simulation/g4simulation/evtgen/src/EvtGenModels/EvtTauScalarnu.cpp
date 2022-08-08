
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

#include "EvtGenModels/EvtTauScalarnu.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <iostream>
#include <stdlib.h>
#include <string>

std::string EvtTauScalarnu::getName()
{
    return "TAUSCALARNU";
}

EvtDecayBase* EvtTauScalarnu::clone()
{
    return new EvtTauScalarnu;
}

void EvtTauScalarnu::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::DIRAC );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::NEUTRINO );
}

void EvtTauScalarnu::initProbMax()
{
    setProbMax( 90.0 );
}

void EvtTauScalarnu::decay( EvtParticle* p )
{
    static EvtId TAUM = EvtPDL::getId( "tau-" );
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle* nut;
    nut = p->getDaug( 1 );
    EvtVector4R momscalar = p->getDaug( 0 )->getP4();

    EvtVector4C tau1, tau2;

    if ( p->getId() == TAUM ) {
        tau1 = EvtLeptonVACurrent( nut->spParentNeutrino(), p->sp( 0 ) );
        tau2 = EvtLeptonVACurrent( nut->spParentNeutrino(), p->sp( 1 ) );
    } else {
        tau1 = EvtLeptonVACurrent( p->sp( 0 ), nut->spParentNeutrino() );
        tau2 = EvtLeptonVACurrent( p->sp( 1 ), nut->spParentNeutrino() );
    }

    vertex( 0, tau1 * momscalar );
    vertex( 1, tau2 * momscalar );

    return;
}
