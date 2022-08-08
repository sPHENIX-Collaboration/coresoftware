
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

#include "EvtGenModels/EvtVSS.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <stdlib.h>
#include <string>

std::string EvtVSS::getName()
{
    return "VSS";
}

EvtDecayBase* EvtVSS::clone()
{
    return new EvtVSS;
}

void EvtVSS::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );

    // check that there are 2 daughters
    checkNDaug( 2 );

    // check the parent and daughter spins
    checkSpinParent( EvtSpinType::VECTOR );
    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
}

void EvtVSS::initProbMax()
{
    setProbMax( 1.0 );
}

void EvtVSS::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R pDaug = p->getDaug( 0 )->getP4();

    double norm = 1.0 / pDaug.d3mag();

    for ( int i = 0; i < 3; i++ )
        vertex( i, norm * pDaug * ( p->eps( i ) ) );

    return;
}
