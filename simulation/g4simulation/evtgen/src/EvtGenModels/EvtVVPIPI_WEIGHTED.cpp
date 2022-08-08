
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

#include "EvtGenModels/EvtVVPIPI_WEIGHTED.hh"

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

std::string EvtVVPIPI_WEIGHTED::getName()
{
    return "VVPIPI_WEIGHTED";
}

EvtDecayBase* EvtVVPIPI_WEIGHTED::clone()
{
    return new EvtVVPIPI_WEIGHTED;
}

void EvtVVPIPI_WEIGHTED::init()
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
            << "EvtVVPIPI_WEIGHTED generator expected "
            << " pi+ and pi- (or pi0 and pi0) "
            << "as 2nd and 3rd daughter. " << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtVVPIPI_WEIGHTED::initProbMax()
{
    //Hard coded... should not be hard to calculate...
    setProbMax( 0.08 * 1.13 );
}

double reweight_event( double pipi_mass )
{
    pipi_mass *= 1000.0;
    return sqrt( -3.6911336508223251 + 0.019119831948029617 * pipi_mass +
                 -1.8962883732377376e-05 * pipi_mass * pipi_mass );
}

void EvtVVPIPI_WEIGHTED::decay( EvtParticle* psi_prime )
{
    psi_prime->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle *jpsi, *pi1, *pi2;

    jpsi = psi_prime->getDaug( 0 );
    pi1 = psi_prime->getDaug( 1 );
    pi2 = psi_prime->getDaug( 2 );

    //  Put phase space results into the daughters.

    EvtVector4C ep0, ep1, ep2;

    ep0 = psi_prime->eps( 0 );
    ep1 = psi_prime->eps( 1 );
    ep2 = psi_prime->eps( 2 );

    EvtVector4C e0, e1, e2;

    e0 = jpsi->epsParent( 0 );
    e1 = jpsi->epsParent( 1 );
    e2 = jpsi->epsParent( 2 );

    double mass2 = ( pi1->getP4() + pi2->getP4() ).mass2();

    double fac = mass2 - 4 * pi1->mass() * pi2->mass();

    fac *= reweight_event( sqrt( mass2 ) );

    vertex( 0, 0, fac * ( ep0 * e0.conj() ) );
    vertex( 0, 1, fac * ( ep0 * e1.conj() ) );
    vertex( 0, 2, fac * ( ep0 * e2.conj() ) );

    vertex( 1, 0, fac * ( ep1 * e0.conj() ) );
    vertex( 1, 1, fac * ( ep1 * e1.conj() ) );
    vertex( 1, 2, fac * ( ep1 * e2.conj() ) );

    vertex( 2, 0, fac * ( ep2 * e0.conj() ) );
    vertex( 2, 1, fac * ( ep2 * e1.conj() ) );
    vertex( 2, 2, fac * ( ep2 * e2.conj() ) );

    return;
}
