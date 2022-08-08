
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

#include "EvtGenModels/EvtSTSCP.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtSTSCP::getName()
{
    return "STS_CP";
}

void EvtSTSCP::initProbMax()
{
    //This is probably not quite right, but it should do as a start...
    //Anders

    setProbMax( 20 * ( getArg( 3 ) * getArg( 3 ) + getArg( 5 ) * getArg( 5 ) ) );
}

EvtDecayBase* EvtSTSCP::clone()
{
    return new EvtSTSCP;
}

void EvtSTSCP::init()
{
    // check that there are 7 arguments
    checkNArg( 7 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::TENSOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
}

void EvtSTSCP::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtParticle* t1;
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    t1 = p->getDaug( 0 );
    EvtVector4R momt = t1->getP4();
    EvtVector4R moms = p->getDaug( 1 )->getP4();
    double masst = t1->mass();

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    double m_parent = p->mass();

    EvtComplex amp;

    EvtComplex A, Abar;

    A = EvtComplex( getArg( 3 ) * cos( getArg( 4 ) ),
                    getArg( 3 ) * sin( getArg( 4 ) ) );
    Abar = EvtComplex( getArg( 5 ) * cos( getArg( 6 ) ),
                       getArg( 5 ) * sin( getArg( 6 ) ) );

    if ( other_b == B0B ) {
        amp = A * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
              EvtComplex( cos( -2.0 * getArg( 0 ) ), sin( -2.0 * getArg( 0 ) ) ) *
                  getArg( 2 ) * EvtComplex( 0.0, 1.0 ) * Abar *
                  sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
    }
    if ( other_b == B0 ) {
        amp = A * EvtComplex( cos( 2.0 * getArg( 0 ) ), sin( 2.0 * getArg( 0 ) ) ) *
                  EvtComplex( 0.0, 1.0 ) *
                  sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) ) +
              getArg( 2 ) * Abar * cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
    }

    EvtVector4R p4_parent;

    p4_parent = momt + moms;

    double norm = masst * masst / ( m_parent * momt.d3mag() * momt.d3mag() );

    vertex( 0, amp * norm * t1->epsTensorParent( 0 ).cont1( p4_parent ) *
                   p4_parent );
    vertex( 1, amp * norm * t1->epsTensorParent( 1 ).cont1( p4_parent ) *
                   p4_parent );
    vertex( 2, amp * norm * t1->epsTensorParent( 2 ).cont1( p4_parent ) *
                   p4_parent );
    vertex( 3, amp * norm * t1->epsTensorParent( 3 ).cont1( p4_parent ) *
                   p4_parent );
    vertex( 4, amp * norm * t1->epsTensorParent( 4 ).cont1( p4_parent ) *
                   p4_parent );

    return;
}

std::string EvtSTSCP::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "weakPhase";
        case 1:
            return "deltaM";
        case 2:
            return "finalStateCP";
        case 3:
            return "Af";
        case 4:
            return "AfPhase";
        case 5:
            return "Abarf";
        case 6:
            return "AbarfPhase";
        default:
            return "";
    }
}

std::string EvtSTSCP::getParamDefault( int i )
{
    switch ( i ) {
        case 3:
            return "1.0";
        case 4:
            return "0.0";
        case 5:
            return "1.0";
        case 6:
            return "0.0";
        default:
            return "";
    }
}
