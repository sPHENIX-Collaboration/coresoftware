
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

#include "EvtGenModels/EvtBTo3piCP.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

std::string EvtBTo3piCP::getName()
{
    return "BTO3PI_CP";
}

EvtBTo3piCP* EvtBTo3piCP::clone()
{
    return new EvtBTo3piCP;
}

void EvtBTo3piCP::init()
{
    // check that there are 2 arguments
    checkNArg( 2 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );
}

void EvtBTo3piCP::initProbMax()
{
    // perform common blocks initialization before
    // first use
    double alpha = getArg( 1 );
    int iset;

    iset = 10000;

    EvtVector4R p4piplus, p4piminus, p4gamm1, p4gamm2;

    double realA, imgA, realbarA, imgbarA;

    generator.Evt3pi( alpha, iset, p4piplus, p4piminus, p4gamm1, p4gamm2, realA,
                      imgA, realbarA, imgbarA );

    setProbMax( 1.5 );
}

void EvtBTo3piCP::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    EvtParticle *pip, *pim, *pi0;

    p->makeDaughters( getNDaug(), getDaugs() );

    //  p->init_daug(SCALAR,&pip,SCALAR,&pim,SCALAR,&pi0);
    pip = p->getDaug( 0 );
    pim = p->getDaug( 1 );
    pi0 = p->getDaug( 2 );

    EvtVector4R p4[3];

    double dm = getArg( 0 );
    double alpha = getArg( 1 );
    int iset;

    iset = 0;

    EvtVector4R p4piplus, p4piminus, p4gamm1, p4gamm2;

    double realA, imgA, realbarA, imgbarA;

    generator.Evt3pi( alpha, iset, p4[0], p4[1], p4gamm1, p4gamm2, realA, imgA,
                      realbarA, imgbarA );

    p4[2] = p4gamm1 + p4gamm2;

    if ( pip->getId() == EvtPDL::getId( "pi+" ) ) {
        pip->init( getDaug( 0 ), p4[0] );
        pim->init( getDaug( 1 ), p4[1] );
    } else {
        pip->init( getDaug( 0 ), p4[1] );
        pim->init( getDaug( 1 ), p4[0] );
    }

    pi0->init( getDaug( 2 ), p4[2] );

    EvtComplex amp;

    EvtComplex A( realA, imgA );
    EvtComplex Abar( realbarA, imgbarA );

    if ( other_b == B0B ) {
        amp = A * cos( dm * t / ( 2 * EvtConst::c ) ) +
              EvtComplex( 0., 1. ) * Abar * sin( dm * t / ( 2 * EvtConst::c ) );
    }
    if ( other_b == B0 ) {
        amp = Abar * cos( dm * t / ( 2 * EvtConst::c ) ) +
              EvtComplex( 0., 1. ) * A * sin( dm * t / ( 2 * EvtConst::c ) );
    }

    vertex( amp );

    return;
}
