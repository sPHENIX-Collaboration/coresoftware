
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

#include "EvtGenModels/EvtSVSNONCPEIGEN.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <stdlib.h>
#include <string>

std::string EvtSVSNONCPEIGEN::getName()
{
    return "SVS_NONCPEIGEN";
}

EvtDecayBase* EvtSVSNONCPEIGEN::clone()
{
    return new EvtSVSNONCPEIGEN;
}

void EvtSVSNONCPEIGEN::init()
{
    // check that there are 11 arguments
    checkNArg( 11, 7 );
    checkNDaug( 2 );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );

    _dm = getArg( 1 );
    _phickm = 2 * getArg( 0 ) + getArg( 2 );

    _A_f = EvtComplex( getArg( 3 ) * cos( getArg( 4 ) ),
                       getArg( 3 ) * sin( getArg( 4 ) ) );
    _Abar_f = EvtComplex( getArg( 5 ) * cos( getArg( 6 ) ),
                          getArg( 5 ) * sin( getArg( 6 ) ) );

    _A_fbar = _Abar_f;
    _Abar_fbar = _A_f;

    if ( getNArg() == 11 ) {
        _A_fbar = EvtComplex( getArg( 7 ) * cos( getArg( 8 ) ),
                              getArg( 7 ) * sin( getArg( 8 ) ) );
        _Abar_fbar = EvtComplex( getArg( 9 ) * cos( getArg( 10 ) ),
                                 getArg( 9 ) * sin( getArg( 10 ) ) );
    }
}

void EvtSVSNONCPEIGEN::initProbMax()
{
    double theProbMax = abs( _A_f ) * abs( _A_f ) +
                        abs( _Abar_f ) * abs( _Abar_f ) +
                        abs( _A_fbar ) * abs( _A_fbar ) +
                        abs( _Abar_fbar ) * abs( _Abar_fbar );

    setProbMax( theProbMax );
}

void EvtSVSNONCPEIGEN::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;
    EvtId daugs[2];

    // MB: flip selects the final of the decay
    int flip = ( ( p->getId() == B0 ) ? 0 : 1 );
    daugs[0] = getDaug( 0 );
    daugs[1] = getDaug( 1 );
    p->initializePhaseSpace( 2, daugs );

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    EvtComplex amp;
    double dmt2 = ( _dm * t ) / ( 2 * EvtConst::c );
    EvtComplex ePlusIPhi( cos( _phickm ), sin( _phickm ) );
    EvtComplex eMinusIPhi( cos( -_phickm ), -sin( _phickm ) );

    // flip == 0 : D-rho+
    // flip == 1 : D+rho-

    if ( !flip ) {
        if ( other_b == B0B ) {
            // At t=0 we have a B0
            amp = cos( dmt2 ) * _A_f +
                  eMinusIPhi * EvtComplex( 0.0, sin( dmt2 ) ) * _Abar_f;
        }
        if ( other_b == B0 ) {
            // At t=0 we have a B0bar
            amp = ePlusIPhi * EvtComplex( 0.0, sin( dmt2 ) ) * _A_f +
                  cos( dmt2 ) * _Abar_f;
        }
    } else {
        if ( other_b == B0B ) {
            // At t=0 we have a B0
            amp = cos( dmt2 ) * _A_fbar +
                  eMinusIPhi * EvtComplex( 0.0, sin( dmt2 ) ) * _Abar_fbar;
        }
        if ( other_b == B0 ) {
            // At t=0 we have a B0bar
            amp = ePlusIPhi * EvtComplex( 0.0, sin( dmt2 ) ) * _A_fbar +
                  cos( dmt2 ) * _Abar_fbar;
        }
    }

    EvtParticle* v;
    v = p->getDaug( 0 );

    EvtVector4R momv = p->getDaug( 0 )->getP4();
    EvtVector4R moms = p->getDaug( 1 )->getP4();
    EvtVector4R p4_parent = momv + moms;

    double norm = momv.mass() / ( momv.d3mag() * p->mass() );

    vertex( 0, amp * norm * p4_parent * ( v->epsParent( 0 ) ) );
    vertex( 1, amp * norm * p4_parent * ( v->epsParent( 1 ) ) );
    vertex( 2, amp * norm * p4_parent * ( v->epsParent( 2 ) ) );

    return;
}
