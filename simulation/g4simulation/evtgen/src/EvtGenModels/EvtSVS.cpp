
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

#include "EvtGenModels/EvtSVS.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtSVS::getName()
{
    return "SVS";
}

EvtDecayBase* EvtSVS::clone()
{
    return new EvtSVS;
}

void EvtSVS::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
}

void EvtSVS::initProbMax()
{
    setProbMax( 1.0 );
}

void EvtSVS::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle* v;
    v = p->getDaug( 0 );
    double massv = v->mass();
    EvtVector4R momv = v->getP4();
    EvtVector4R moms = p->getDaug( 1 )->getP4();
    double m_parent = p->mass();
    EvtVector4R p4_parent = momv + moms;

    double norm = massv / ( momv.d3mag() * m_parent );
    p4_parent = norm * p4_parent;
    vertex( 0, p4_parent * ( v->epsParent( 0 ) ) );
    vertex( 1, p4_parent * ( v->epsParent( 1 ) ) );
    vertex( 2, p4_parent * ( v->epsParent( 2 ) ) );

    return;
}
