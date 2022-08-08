
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

#include "EvtGenModels/EvtVSPPwave.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtVSPPwave::getName()
{
    return "VSP_PWAVE";
}

EvtDecayBase* EvtVSPPwave::clone()
{
    return new EvtVSPPwave;
}

void EvtVSPPwave::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::VECTOR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::PHOTON );
}

void EvtVSPPwave::initProbMax()
{
    setProbMax( 1 );
}

void EvtVSPPwave::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle* gamma;
    gamma = p->getDaug( 1 );

    double m_p = p->mass();
    EvtVector4R momgamma = gamma->getP4();

    //work in the parent ,p,  rest frame.
    EvtVector4R p4_p;
    p4_p.set( m_p, 0.0, 0.0, 0.0 );

    //  Put phase space results into the daughters.

    EvtTensor4C tds;

    double norm = 1 / ( m_p * momgamma.d3mag() );

    tds = dual( EvtGenFunctions::directProd( norm * p4_p, momgamma ) );

    vertex(
        0, 0,
        ( tds.cont1( p->eps( 0 ) ) ).cont( gamma->epsParentPhoton( 0 ).conj() ) );
    vertex(
        0, 1,
        ( tds.cont1( p->eps( 0 ) ) ).cont( gamma->epsParentPhoton( 1 ).conj() ) );

    vertex(
        1, 0,
        ( tds.cont1( p->eps( 1 ) ) ).cont( gamma->epsParentPhoton( 0 ).conj() ) );
    vertex(
        1, 1,
        ( tds.cont1( p->eps( 1 ) ) ).cont( gamma->epsParentPhoton( 1 ).conj() ) );

    vertex(
        2, 0,
        ( tds.cont1( p->eps( 2 ) ) ).cont( gamma->epsParentPhoton( 0 ).conj() ) );
    vertex(
        2, 1,
        ( tds.cont1( p->eps( 2 ) ) ).cont( gamma->epsParentPhoton( 1 ).conj() ) );

    return;
}
