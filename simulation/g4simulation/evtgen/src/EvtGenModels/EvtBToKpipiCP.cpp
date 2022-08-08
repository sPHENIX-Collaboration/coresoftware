
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

#include "EvtGenModels/EvtBToKpipiCP.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <stdlib.h>
#include <string>

std::string EvtBToKpipiCP::getName()
{
    return "BTOKPIPI_CP";
}

EvtBToKpipiCP* EvtBToKpipiCP::clone()
{
    return new EvtBToKpipiCP;
}

void EvtBToKpipiCP::init()
{
    // check that there are 3 arguments
    checkNArg( 3 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );

    double alpha = getArg( 1 );
    double beta = getArg( 2 );
    int iset;
    iset = 10000;

    EvtVector4R p4Kplus, p4piminus, p4gamm1, p4gamm2;

    double realA, imgA, realbarA, imgbarA;

    generator.EvtKpipi( alpha, beta, iset, p4Kplus, p4piminus, p4gamm1, p4gamm2,
                        realA, imgA, realbarA, imgbarA );
}

void EvtBToKpipiCP::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    EvtParticle *Kp, *pim, *pi0;

    p->makeDaughters( getNDaug(), getDaugs() );
    Kp = p->getDaug( 0 );
    pim = p->getDaug( 1 );
    pi0 = p->getDaug( 2 );

    EvtVector4R p4[3];

    //double dm=getArg(0);
    double alpha = getArg( 1 );
    double beta = getArg( 2 );
    int iset;

    iset = 0;

    EvtVector4R p4Kplus, p4piminus, p4gamm1, p4gamm2;

    double realA, imgA, realbarA, imgbarA;

    generator.EvtKpipi( alpha, beta, iset, p4[0], p4[1], p4gamm1, p4gamm2,
                        realA, imgA, realbarA, imgbarA );

    p4[2] = p4gamm1 + p4gamm2;

    Kp->init( getDaug( 0 ), p4[0] );
    pim->init( getDaug( 1 ), p4[1] );
    pi0->init( getDaug( 2 ), p4[2] );

    EvtComplex amp;

    EvtComplex A( realA, imgA );
    EvtComplex Abar( realbarA, imgbarA );

    if ( other_b == B0B ) {
        amp = Abar;
    }
    if ( other_b == B0 ) {
        amp = A;
    }

    vertex( amp );

    return;
}
