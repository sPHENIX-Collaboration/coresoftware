
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

#include "EvtGenModels/EvtSTS.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtSTS::getName()
{
    return "STS";
}

EvtDecayBase* EvtSTS::clone()
{
    return new EvtSTS;
}

void EvtSTS::initProbMax()
{
    setProbMax( 20.0 );
}

void EvtSTS::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::TENSOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
}

void EvtSTS::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle* t1 = p->getDaug( 0 );

    EvtVector4R momt = t1->getP4();
    EvtVector4R moms = p->getDaug( 1 )->getP4();
    double masst = t1->mass();

    EvtVector4R p4_parent = momt + moms;

    double m_parent = p4_parent.mass();

    double norm = masst * masst / ( m_parent * momt.d3mag() * momt.d3mag() );

    vertex( 0, norm * t1->epsTensorParent( 0 ).cont1( p4_parent ) * p4_parent );
    vertex( 1, norm * t1->epsTensorParent( 1 ).cont1( p4_parent ) * p4_parent );
    vertex( 2, norm * t1->epsTensorParent( 2 ).cont1( p4_parent ) * p4_parent );
    vertex( 3, norm * t1->epsTensorParent( 3 ).cont1( p4_parent ) * p4_parent );
    vertex( 4, norm * t1->epsTensorParent( 4 ).cont1( p4_parent ) * p4_parent );

    return;
}
