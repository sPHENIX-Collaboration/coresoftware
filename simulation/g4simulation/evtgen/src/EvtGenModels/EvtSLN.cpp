
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

#include "EvtGenModels/EvtSLN.hh"

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

std::string EvtSLN::getName()
{
    return "SLN";
}

EvtDecayBase* EvtSLN::clone()
{
    return new EvtSLN;
}

void EvtSLN::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::NEUTRINO );
}

void EvtSLN::initProbMax()
{
    double M = EvtPDL::getMeanMass( getParentId() );
    double m = EvtPDL::getMeanMass( getDaug( 0 ) );

    double probMax = 8.0 * ( M * M - m * m ) * m * m;

    setProbMax( probMax );
}

void EvtSLN::decay( EvtParticle* p )
{
    static EvtId EM = EvtPDL::getId( "e-" );
    static EvtId MUM = EvtPDL::getId( "mu-" );
    static EvtId TAUM = EvtPDL::getId( "tau-" );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle *l, *nul;
    l = p->getDaug( 0 );
    nul = p->getDaug( 1 );

    EvtVector4R p4_p;
    p4_p.set( p->mass(), 0.0, 0.0, 0.0 );

    EvtVector4C l1, l2;

    if ( getDaug( 0 ) == TAUM || getDaug( 0 ) == MUM || getDaug( 0 ) == EM ) {
        l1 = EvtLeptonVACurrent( l->spParent( 0 ), nul->spParentNeutrino() );
        l2 = EvtLeptonVACurrent( l->spParent( 1 ), nul->spParentNeutrino() );
    } else {
        l1 = EvtLeptonVACurrent( nul->spParentNeutrino(), l->spParent( 0 ) );
        l2 = EvtLeptonVACurrent( nul->spParentNeutrino(), l->spParent( 1 ) );
    }

    vertex( 0, p4_p * l1 );
    vertex( 1, p4_p * l2 );

    return;
}
