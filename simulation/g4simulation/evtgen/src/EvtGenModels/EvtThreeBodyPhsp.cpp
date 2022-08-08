
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

#include "EvtGenModels/EvtThreeBodyPhsp.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

std::string EvtThreeBodyPhsp::getName()
{
    return "THREEBODYPHSP";
}

EvtDecayBase* EvtThreeBodyPhsp::clone()
{
    return new EvtThreeBodyPhsp;
}

void EvtThreeBodyPhsp::init()
{
    // check correct number of arguments
    checkNArg( 2, 4 );
    // check number of daughters
    checkNDaug( 3 );
    // Get box size
    m_m12SqMin = getArg( 0 );
    m_m12SqMax = getArg( 1 );
    if ( getNArg() > 2 ) {
        m_m23SqMin = getArg( 2 );
        m_m23SqMax = getArg( 3 );
    } else {
        m_m23SqMin = 0;
        m_m23SqMax = 1e10;
    }
}

void EvtThreeBodyPhsp::initProbMax()
{
    noProbMax();
}

void EvtThreeBodyPhsp::decay( EvtParticle* p )
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

    const double m12borderMin = mDaug1 + mDaug2;
    const double m12borderMax = mParent - mDaug3;
    const double m12Min = std::max( m_m12SqMin, m12borderMin * m12borderMin );
    const double m12Max = std::min( m_m12SqMax, m12borderMax * m12borderMax );

    const double m23borderMin = mDaug2 + mDaug3;
    const double m23borderMax = mParent - mDaug1;
    const double m23Min = std::max( m_m23SqMin, m23borderMin * m23borderMin );
    const double m23Max = std::min( m_m23SqMax, m23borderMax * m23borderMax );

    const int nMaxTrials = 1000;
    int iTrial = 0;
    bool goodEvent = false;
    double m12Sq, m23Sq, m23LowLimit, m23HighLimit;
    do {
        m12Sq = EvtRandom::Flat( m12Min, m12Max );
        m23Sq = EvtRandom::Flat( m23Min, m23Max );
        // By definition, m12Sq is always within Dalitz plot, but need to check
        // that also m23Sq is in
        double E2st = 0.5 * ( m12Sq - mDaug1 * mDaug1 + mDaug2 * mDaug2 ) /
                      std::sqrt( m12Sq );
        double E3st = 0.5 * ( mParent * mParent - m12Sq - mDaug3 * mDaug3 ) /
                      std::sqrt( m12Sq );
        double p2st = std::sqrt( E2st * E2st - mDaug2 * mDaug2 );
        double p3st = std::sqrt( E3st * E3st - mDaug3 * mDaug3 );
        m23HighLimit = ( E2st + E3st ) * ( E2st + E3st ) -
                       ( p2st - p3st ) * ( p2st - p3st );
        m23LowLimit = ( E2st + E3st ) * ( E2st + E3st ) -
                      ( p2st + p3st ) * ( p2st + p3st );
        if ( m23Sq > m23LowLimit && m23Sq < m23HighLimit ) {
            goodEvent = true;
        }
        ++iTrial;
    } while ( goodEvent == false && iTrial < nMaxTrials );
    if ( !goodEvent ) {
        EvtGenReport( EVTGEN_WARNING, "EvtThreeBodyPhsp" )
            << "Failed to generate m12Sq and m23Sq. Taking last m12Sq and midpoint of allowed m23Sq.\n";
        m23Sq = 0.5 * ( m23LowLimit + m23HighLimit );
    }

    // At this moment we have valid invariant masses squared
    EvtGenKine::ThreeBodyKine( m12Sq, m23Sq, p );

    return;
}
