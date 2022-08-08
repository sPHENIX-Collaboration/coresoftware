
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

#include "EvtGenModels/EvtVVpipi.hh"

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

std::string EvtVVpipi::getName()
{
    return "VVPIPI";
}

EvtDecayBase* EvtVVpipi::clone()
{
    return new EvtVVpipi;
}

void EvtVVpipi::init()
{
    static EvtId PIP = EvtPDL::getId( "pi+" );
    static EvtId PIM = EvtPDL::getId( "pi-" );
    static EvtId PI0 = EvtPDL::getId( "pi0" );

    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::VECTOR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );

    if ( ( !( getDaug( 1 ) == PIP && getDaug( 2 ) == PIM ) ) &&
         ( !( getDaug( 1 ) == PI0 && getDaug( 2 ) == PI0 ) ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtVVpipi generator expected "
            << " pi+ and pi- (or pi0 and pi0) "
            << "as 2nd and 3rd daughter. " << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtVVpipi::initProbMax()
{
    //Hard coded... should not be hard to calculate...
    setProbMax( 0.08 );
}

void EvtVVpipi::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle *v, *s1, *s2;

    v = p->getDaug( 0 );
    s1 = p->getDaug( 1 );
    s2 = p->getDaug( 2 );

    //  Put phase space results into the daughters.

    EvtVector4C ep0, ep1, ep2;

    ep0 = p->eps( 0 );
    ep1 = p->eps( 1 );
    ep2 = p->eps( 2 );

    double fac = ( s1->getP4() + s2->getP4() ).mass2() -
                 4 * s1->mass() * s2->mass();

    vertex( 0, 0, fac * ( ep0 * v->epsParent( 0 ).conj() ) );
    vertex( 0, 1, fac * ( ep0 * v->epsParent( 1 ).conj() ) );
    vertex( 0, 2, fac * ( ep0 * v->epsParent( 2 ).conj() ) );

    vertex( 1, 0, fac * ( ep1 * v->epsParent( 0 ).conj() ) );
    vertex( 1, 1, fac * ( ep1 * v->epsParent( 1 ).conj() ) );
    vertex( 1, 2, fac * ( ep1 * v->epsParent( 2 ).conj() ) );

    vertex( 2, 0, fac * ( ep2 * v->epsParent( 0 ).conj() ) );
    vertex( 2, 1, fac * ( ep2 * v->epsParent( 1 ).conj() ) );
    vertex( 2, 2, fac * ( ep2 * v->epsParent( 2 ).conj() ) );

    return;
}
