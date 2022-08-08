
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

#include "EvtGenModels/EvtOmegaDalitz.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtOmegaDalitz::getName()
{
    return "OMEGA_DALITZ";
}

EvtDecayBase* EvtOmegaDalitz::clone()
{
    return new EvtOmegaDalitz;
}

void EvtOmegaDalitz::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::VECTOR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );
}

void EvtOmegaDalitz::initProbMax()
{
    setProbMax( 1.0 );
}

void EvtOmegaDalitz::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4C ep[3];

    ep[0] = p->eps( 0 );
    ep[1] = p->eps( 1 );
    ep[2] = p->eps( 2 );

    EvtVector4R mompi1 = p->getDaug( 0 )->getP4();
    EvtVector4R mompi2 = p->getDaug( 1 )->getP4();

    EvtVector3R p1( mompi1.get( 1 ), mompi1.get( 2 ), mompi1.get( 3 ) );
    EvtVector3R p2( mompi2.get( 1 ), mompi2.get( 2 ), mompi2.get( 3 ) );
    EvtVector3R q = cross( p2, p1 );

    EvtVector3C e1( ep[0].get( 1 ), ep[0].get( 2 ), ep[0].get( 3 ) );
    EvtVector3C e2( ep[1].get( 1 ), ep[1].get( 2 ), ep[1].get( 3 ) );
    EvtVector3C e3( ep[2].get( 1 ), ep[2].get( 2 ), ep[2].get( 3 ) );

    //This is an approximate formula of the maximum value that
    //|q| can have.
    double norm = 1.14 / ( p->mass() * p->mass() / 9.0 - mompi1.mass2() );

    vertex( 0, norm * e1 * q );
    vertex( 1, norm * e2 * q );
    vertex( 2, norm * e3 * q );

    return;
}
