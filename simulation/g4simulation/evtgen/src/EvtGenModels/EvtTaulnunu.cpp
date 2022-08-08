
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

#include "EvtGenModels/EvtTaulnunu.hh"

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

std::string EvtTaulnunu::getName()
{
    return "TAULNUNU";
}

EvtDecayBase* EvtTaulnunu::clone()
{
    return new EvtTaulnunu;
}

void EvtTaulnunu::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::DIRAC );

    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::NEUTRINO );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );
}

void EvtTaulnunu::initProbMax()
{
    setProbMax( 650.0 );
}

void EvtTaulnunu::decay( EvtParticle* p )
{
    static EvtId TAUM = EvtPDL::getId( "tau-" );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle *l, *nul, *nut;

    l = p->getDaug( 0 );
    nul = p->getDaug( 1 );
    nut = p->getDaug( 2 );

    EvtVector4C l1, l2, tau1, tau2;

    if ( p->getId() == TAUM ) {
        tau1 = EvtLeptonVACurrent( nut->spParentNeutrino(), p->sp( 0 ) );
        tau2 = EvtLeptonVACurrent( nut->spParentNeutrino(), p->sp( 1 ) );
        l1 = EvtLeptonVACurrent( l->spParent( 0 ), nul->spParentNeutrino() );
        l2 = EvtLeptonVACurrent( l->spParent( 1 ), nul->spParentNeutrino() );

    } else {
        tau1 = EvtLeptonVACurrent( p->sp( 0 ), nut->spParentNeutrino() );
        tau2 = EvtLeptonVACurrent( p->sp( 1 ), nut->spParentNeutrino() );
        l1 = EvtLeptonVACurrent( nul->spParentNeutrino(), l->spParent( 0 ) );
        l2 = EvtLeptonVACurrent( nul->spParentNeutrino(), l->spParent( 1 ) );
    }

    vertex( 0, 0, tau1 * l1 );
    vertex( 0, 1, tau1 * l2 );
    vertex( 1, 0, tau2 * l1 );
    vertex( 1, 1, tau2 * l2 );

    return;
}
