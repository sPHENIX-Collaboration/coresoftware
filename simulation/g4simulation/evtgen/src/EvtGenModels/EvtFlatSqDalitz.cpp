
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

#include "EvtGenModels/EvtFlatSqDalitz.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <cmath>
#include <string>

EvtFlatSqDalitz::~EvtFlatSqDalitz()
{
}

std::string EvtFlatSqDalitz::getName()
{
    return "FLATSQDALITZ";
}

EvtDecayBase* EvtFlatSqDalitz::clone()
{
    return new EvtFlatSqDalitz;
}

void EvtFlatSqDalitz::initProbMax()
{
    noProbMax();
}

void EvtFlatSqDalitz::init()
{
    //check there are 3 daughters
    checkNDaug( 3 );

    // check that there are 0 arguments
    checkNArg( 0, 2, 4 );

    if ( getNArg() > 0 ) {
        m_mPrimeMin = getArg( 0 );
        m_mPrimeMax = getArg( 1 );
    }
    if ( getNArg() > 2 ) {
        m_thetaPrimeMin = getArg( 2 );
        m_thetaPrimeMax = getArg( 3 );
    }
}

void EvtFlatSqDalitz::decay( EvtParticle* p )
{
    p->makeDaughters( getNDaug(), getDaugs() );
    p->generateMassTree();
    const double mParent = p->mass();
    EvtParticle* daug1 = p->getDaug( 0 );
    EvtParticle* daug2 = p->getDaug( 1 );
    EvtParticle* daug3 = p->getDaug( 2 );
    const double mDaug1 = daug1->mass();
    const double mDaug2 = daug2->mass();
    const double mDaug3 = daug3->mass();
    const double mParentSq = mParent * mParent;
    const double mDaug1Sq = mDaug1 * mDaug1;
    const double mDaug2Sq = mDaug2 * mDaug2;
    const double mDaug3Sq = mDaug3 * mDaug3;

    // Generate m' and theta'
    const double mPrime = EvtRandom::Flat( m_mPrimeMin, m_mPrimeMax );
    const double thetaPrime = EvtRandom::Flat( m_thetaPrimeMin, m_thetaPrimeMax );

    // calculate m12 and m23
    const double m12 = 0.5 * ( std::cos( mPrime * EvtConst::pi ) + 1 ) *
                           ( mParent - ( mDaug1 + mDaug2 + mDaug3 ) ) +
                       mDaug1 + mDaug2;
    const double m12Sq = m12 * m12;

    const double en1 = ( m12Sq - mDaug2Sq + mDaug1Sq ) / ( 2. * m12 );
    const double en3 = ( mParentSq - m12Sq - mDaug3Sq ) / ( 2. * m12 );

    const double p1 = std::sqrt( en1 * en1 - mDaug1Sq );
    const double p3 = std::sqrt( en3 * en3 - mDaug3Sq );
    const double m13Sq =
        mDaug1Sq + mDaug3Sq +
        2.0 * ( en1 * en3 - p1 * p3 * std::cos( EvtConst::pi * thetaPrime ) );
    const double m23Sq = mParentSq - m12Sq - m13Sq + mDaug1Sq + mDaug2Sq +
                         mDaug3Sq;

    // Turn m12 and m23 into momenta
    EvtGenKine::ThreeBodyKine( m12Sq, m23Sq, p );

    return;
}
