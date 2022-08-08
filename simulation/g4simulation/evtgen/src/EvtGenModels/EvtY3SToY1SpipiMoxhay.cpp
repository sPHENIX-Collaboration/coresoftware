
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

#include "EvtGenModels/EvtY3SToY1SpipiMoxhay.hh"

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

std::string EvtY3SToY1SpipiMoxhay::getName()
{
    return "Y3STOY1SPIPIMOXHAY";
}

EvtDecayBase* EvtY3SToY1SpipiMoxhay::clone()
{
    return new EvtY3SToY1SpipiMoxhay;
}

void EvtY3SToY1SpipiMoxhay::init()
{
    static EvtId PIP = EvtPDL::getId( "pi+" );
    static EvtId PIM = EvtPDL::getId( "pi-" );
    static EvtId PI0 = EvtPDL::getId( "pi0" );

    // check that there are 2 arguments
    checkNArg( 2 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::VECTOR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );

    if ( ( !( getDaug( 1 ) == PIP && getDaug( 2 ) == PIM ) ) &&
         ( !( getDaug( 1 ) == PI0 && getDaug( 2 ) == PI0 ) ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtY3SToY1SpipiMoxhay generator expected "
            << " pi+ and pi- (or pi0 and pi0) "
            << "as 2nd and 3rd daughter. " << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtY3SToY1SpipiMoxhay::initProbMax()
{
    setProbMax( 0.2 );
}

void EvtY3SToY1SpipiMoxhay::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtParticle *v, *s1, *s2;

    v = p->getDaug( 0 );
    s1 = p->getDaug( 1 );
    s2 = p->getDaug( 2 );

    // setup the parameters needed for this model
    double g_spp = 0.64;
    double lambda = -0.73;
    double m_sigma = 0.71;
    double f_pi = 0.094;
    double m_pi = s1->getP4().mass();

    double MV1 = p->getP4().mass();
    double MV2 = v->getP4().mass();

    double q = ( s1->getP4() + s2->getP4() ).mass();

    double EV2 = ( MV1 * MV1 - MV2 * MV2 - q * q ) / ( 2.0 * q );

    double ReB_over_A = getArg( 0 );
    double ImB_over_A = getArg( 1 );

    EvtComplex Xi;

    Xi = EvtComplex( 2.0 / EvtConst::pi *
                         ( 1.0 - sqrt( 1.0 - 4 * m_pi * m_pi / ( q * q ) ) *
                                     log( ( sqrt( q * q ) +
                                            sqrt( q * q - 4.0 * m_pi * m_pi ) ) /
                                          ( 2 * m_pi ) ) ),
                     sqrt( 1.0 - 4 * m_pi * m_pi / ( q * q ) ) );

    // The form factor
    EvtComplex F;

    F = ( g_spp * g_spp + lambda * ( m_sigma * m_sigma - q * q ) ) /
        ( ( ( m_sigma * m_sigma - q * q ) * ( 1.0 - lambda * Xi ) -
            ( g_spp * g_spp * Xi ) ) *
          1.0 / ( 8.0 * EvtConst::pi * f_pi * f_pi ) * q * q );

    EvtComplex B_over_A;
    B_over_A = EvtComplex( ReB_over_A, ImB_over_A );

    // The dGamma/d(M_pipi) spectrum
    EvtComplex dGdMpp;

    dGdMpp = abs2( ( q * q * F - B_over_A ) ) * q *
             sqrt( q * q - 4 * m_pi * m_pi ) * sqrt( EV2 * EV2 - MV2 * MV2 );

    setProb( real( dGdMpp ) );
    return;
}
