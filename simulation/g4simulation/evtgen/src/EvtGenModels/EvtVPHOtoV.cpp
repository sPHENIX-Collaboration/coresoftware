
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

#include "EvtGenModels/EvtVPHOtoV.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <stdlib.h>
#include <string>

std::string EvtVPHOtoV::getName()
{
    return "VPHOTOV";
}

EvtDecayBase* EvtVPHOtoV::clone()
{
    return new EvtVPHOtoV;
}

void EvtVPHOtoV::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );

    // check that there are 1 daughters
    checkNDaug( 1 );

    // check the parent and daughter spins
    checkSpinParent( EvtSpinType::VECTOR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );
}

void EvtVPHOtoV::initProbMax()
{
    setProbMax( 1.0 );
}

void EvtVPHOtoV::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle* d = p->getDaug( 0 );

    d->setP4( p->getP4Restframe() );

    vertex( 0, 0, p->eps( 0 ) * p->epsParent( 0 ).conj() );
    vertex( 1, 0, p->eps( 1 ) * p->epsParent( 0 ).conj() );
    vertex( 2, 0, p->eps( 2 ) * p->epsParent( 0 ).conj() );

    vertex( 0, 1, p->eps( 0 ) * p->epsParent( 1 ).conj() );
    vertex( 1, 1, p->eps( 1 ) * p->epsParent( 1 ).conj() );
    vertex( 2, 1, p->eps( 2 ) * p->epsParent( 1 ).conj() );

    vertex( 0, 2, p->eps( 0 ) * p->epsParent( 2 ).conj() );
    vertex( 1, 2, p->eps( 1 ) * p->epsParent( 2 ).conj() );
    vertex( 2, 2, p->eps( 2 ) * p->epsParent( 2 ).conj() );

    return;
}
