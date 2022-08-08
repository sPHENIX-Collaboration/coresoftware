
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

#include "EvtGenModels/EvtTSS.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtTSS::getName()
{
    return "TSS";
}

EvtDecayBase* EvtTSS::clone()
{
    return new EvtTSS;
}

void EvtTSS::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );

    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::TENSOR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
}

void EvtTSS::initProbMax()
{
    setProbMax( 1.0 );
}

void EvtTSS::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R moms1 = p->getDaug( 0 )->getP4();

    double norm = 1.0 / ( moms1.d3mag() * moms1.d3mag() );

    vertex( 0, norm * ( p->epsTensor( 0 ).cont1( EvtVector4C( moms1 ) ) *
                        ( moms1 ) ) );
    vertex( 1, norm * ( p->epsTensor( 1 ).cont1( EvtVector4C( moms1 ) ) *
                        ( moms1 ) ) );
    vertex( 2, norm * ( p->epsTensor( 2 ).cont1( EvtVector4C( moms1 ) ) *
                        ( moms1 ) ) );
    vertex( 3, norm * ( p->epsTensor( 3 ).cont1( EvtVector4C( moms1 ) ) *
                        ( moms1 ) ) );
    vertex( 4, norm * ( p->epsTensor( 4 ).cont1( EvtVector4C( moms1 ) ) *
                        ( moms1 ) ) );

    return;
}
