
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

#include "EvtGenModels/EvtCBTo3piMPP.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

std::string EvtCBTo3piMPP::getName()
{
    return "CB3PI-MPP";
}

EvtCBTo3piMPP* EvtCBTo3piMPP::clone()
{
    return new EvtCBTo3piMPP;
}

void EvtCBTo3piMPP::init()
{
    // check that there are 1 argument
    checkNArg( 1 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );

    EvtVector4R p4[3];
    double alpha = getArg( 0 );

    int iset( 10000 );
    double realA, imgA, realbarA, imgbarA;

    generator.Evt3piMPP( alpha, iset, p4[0], p4[1], p4[2], realA, imgA,
                         realbarA, imgbarA );
}

void EvtCBTo3piMPP::initProbMax()
{
    setProbMax( 1.5 );
}

void EvtCBTo3piMPP::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId BM = EvtPDL::getId( "B-" );
    static EvtId BP = EvtPDL::getId( "B+" );

    EvtParticle *pi1, *pi2, *pi3;

    p->makeDaughters( getNDaug(), getDaugs() );
    pi1 = p->getDaug( 0 );
    pi2 = p->getDaug( 1 );
    pi3 = p->getDaug( 2 );

    EvtVector4R p4[3];
    double alpha = getArg( 0 );

    int iset( 0 );
    double realA, imgA, realbarA, imgbarA;

    generator.Evt3piMPP( alpha, iset, p4[0], p4[1], p4[2], realA, imgA,
                         realbarA, imgbarA );

    pi1->init( getDaug( 0 ), p4[0] );
    pi2->init( getDaug( 1 ), p4[1] );
    pi3->init( getDaug( 2 ), p4[2] );

    EvtComplex A( realA, imgA );
    EvtComplex Abar( realbarA, imgbarA );

    //amp is filled just to make sure the compiler will
    //do its job!! but one has to define amp differently
    // if one wants the B+ or the B- to decay to 3pi!
    //

    EvtComplex amp;
    if ( p->getId() == BP ) {
        amp = A;
    }
    if ( p->getId() == BM ) {
        amp = Abar;
    }

    vertex( amp );

    return;
}
