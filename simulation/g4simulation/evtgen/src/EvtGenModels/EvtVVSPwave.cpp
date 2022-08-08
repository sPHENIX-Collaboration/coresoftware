
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

#include "EvtGenModels/EvtVVSPwave.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtVVSPwave::getName()
{
    return "VVS_PWAVE";
}

EvtDecayBase* EvtVVSPwave::clone()
{
    return new EvtVVSPwave;
}

void EvtVVSPwave::init()
{
    // check that there are 6 arguments
    checkNArg( 6 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::VECTOR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
}

void EvtVVSPwave::initProbMax()
{
    //probmax is 1.0 for all possible decays I think!

    setProbMax( 1.0 );
}

void EvtVVSPwave::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtComplex as( getArg( 0 ) * cos( getArg( 1 ) ),
                   getArg( 0 ) * sin( getArg( 1 ) ) );
    EvtComplex ap( getArg( 2 ) * cos( getArg( 3 ) ),
                   getArg( 2 ) * sin( getArg( 3 ) ) );
    EvtComplex ad( getArg( 4 ) * cos( getArg( 5 ) ),
                   getArg( 4 ) * sin( getArg( 5 ) ) );

    if ( ap != EvtComplex( 0.0, 0.0 ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "In EvtVectorToVectorScalar.cc" << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "P wave not yet implemented!!" << endl;
        ::abort();
    }

    EvtParticle* v;
    v = p->getDaug( 0 );

    EvtTensor4C d, g;

    g.setdiag( 1.0, -1.0, -1.0, -1.0 );

    d = ad * ( ( 1.0 / ( v->getP4().d3mag() * v->getP4().d3mag() ) ) *
                   EvtGenFunctions::directProd( v->getP4(), v->getP4() ) +
               ( 1 / 3.0 ) * g ) +
        as * g;

    EvtVector4C ep0, ep1, ep2;

    ep0 = d.cont1( p->eps( 0 ) );
    ep1 = d.cont1( p->eps( 1 ) );
    ep2 = d.cont1( p->eps( 2 ) );

    vertex( 0, 0, ep0.cont( v->eps( 0 ).conj() ) );
    vertex( 0, 1, ep0.cont( v->eps( 1 ).conj() ) );
    vertex( 0, 2, ep0.cont( v->eps( 2 ).conj() ) );

    vertex( 1, 0, ep1.cont( v->eps( 0 ).conj() ) );
    vertex( 1, 1, ep1.cont( v->eps( 1 ).conj() ) );
    vertex( 1, 2, ep1.cont( v->eps( 2 ).conj() ) );

    vertex( 2, 0, ep2.cont( v->eps( 0 ).conj() ) );
    vertex( 2, 1, ep2.cont( v->eps( 1 ).conj() ) );
    vertex( 2, 2, ep2.cont( v->eps( 2 ).conj() ) );

    return;
}
