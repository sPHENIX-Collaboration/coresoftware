
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

#include "EvtGenModels/EvtPi0Dalitz.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using std::fstream;

std::string EvtPi0Dalitz::getName()
{
    return "PI0_DALITZ";
}

EvtDecayBase* EvtPi0Dalitz::clone()
{
    return new EvtPi0Dalitz;
}

void EvtPi0Dalitz::initProbMax()
{
    // Search for maximum probability. In order to avoid setting up all
    // particles with spinors and four-momenta, we use result after all
    // contractions, which is:
    // 1/((m_R^2-q^2)^2+m_R^2 Gamma_R^2) 1/(2q^2) (M^2-q^2)^2 beta_l
    // (1+cos(theta)^2) where we set cos(theta)=1
    auto daughter1 = getDaug( 0 );
    auto daughter2 = getDaug( 1 );
    double q2Min = EvtPDL::getMass( daughter1 ) + EvtPDL::getMass( daughter2 );
    q2Min *= q2Min;
    double q2Max = EvtPDL::getMass( getParentId() );
    q2Max *= q2Max;
    const int steps = 20000;
    const double step = ( q2Max - q2Min ) / steps;
    double maxProb = 0;
    for ( int ii = 0; ii < steps; ++ii ) {
        double q2 = q2Min + ii * step;
        const double mSqDiff = m_m0Sq - q2;
        const double q2Sq = q2 * q2;
        double prob = ( q2Max - q2 ) * ( q2Max - q2 ) * ( q2 - q2Min ) /
                      ( q2Sq );
        prob *= ( 1.0 / ( mSqDiff * mSqDiff + m_m0SqG0Sq ) );
        // When generating events, we do not start from phase-space, but
        // add some pole to it, weight of which is taken into account
        // elsewhere
        prob /= 1.0 + m_poleSize / ( q2Sq );
        if ( prob > maxProb ) {
            maxProb = prob;
        }
    }
    setProbMax( maxProb * 1.05  );
}

void EvtPi0Dalitz::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
    checkNDaug( 3 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::PHOTON );

    // Rescale pole size to improve efficiency.  Not sure about exact
    // factor, but this seem to be best simple rescaling for
    // eta-->e+e-gamma.
    const double parentMass = EvtPDL::getMass( getParentId() );
    m_poleSize *= parentMass * parentMass / ( 0.135 * 0.135 );
}

void EvtPi0Dalitz::decay( EvtParticle* p )
{
    EvtParticle *ep, *em, *gamma;
    setWeight( p->initializePhaseSpace( getNDaug(), getDaugs(), false,
                                        m_poleSize, 0, 1 ) );
    ep = p->getDaug( 0 );
    em = p->getDaug( 1 );
    gamma = p->getDaug( 2 );

    // the next four lines generates events with a weight such that
    // the efficiency for selecting them is good. The parameter below of
    // 0.1 is the size of the peak at low q^2 (in arbitrary units).
    // The value of 0.1 is appropriate for muons.
    // when you use this remember to remove the cut on q^2!

    //ep em invariant mass^2
    double m2 = ( ep->getP4() + em->getP4() ).mass2();
    EvtVector4R q = ep->getP4() + em->getP4();
    //Just use the prob summed over spins...

    EvtTensor4C w, v;

    v = 2.0 * ( gamma->getP4() * q ) *
            EvtGenFunctions::directProd( q, gamma->getP4() ) -
        ( gamma->getP4() * q ) * ( gamma->getP4() * q ) * EvtTensor4C::g() -
        m2 * EvtGenFunctions::directProd( gamma->getP4(), gamma->getP4() );

    w = 4.0 * ( EvtGenFunctions::directProd( ep->getP4(), em->getP4() ) +
                EvtGenFunctions::directProd( em->getP4(), ep->getP4() ) -
                EvtTensor4C::g() *
                    ( ep->getP4() * em->getP4() - ep->getP4().mass2() ) );

    double prob = ( real( cont( v, w ) ) ) / ( m2 * m2 );
    const double m2Diff = m_m0Sq - m2;
    prob *= ( 1.0 / ( m2Diff * m2Diff + m_m0SqG0Sq ) );

    //  EvtGenReport(EVTGEN_INFO,"EvtGen") << "prob is "<<prob<<endl;
    setProb( prob );

    return;
}
