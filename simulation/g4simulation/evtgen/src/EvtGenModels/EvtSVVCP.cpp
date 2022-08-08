
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

#include "EvtGenModels/EvtSVVCP.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtSVVHelAmp.hh"

#include <stdlib.h>
#include <string>

std::string EvtSVVCP::getName()
{
    return "SVV_CP";
}

EvtDecayBase* EvtSVVCP::clone()
{
    return new EvtSVVCP;
}

void EvtSVVCP::init()
{
    // check that there are 9 arguments
    checkNArg( 9 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::VECTOR );
}

void EvtSVVCP::initProbMax()
{
    //This is probably not quite right, but it should do as a start...
    //Anders

    setProbMax( 2 * ( getArg( 3 ) * getArg( 3 ) + getArg( 5 ) * getArg( 5 ) +
                      getArg( 7 ) * getArg( 7 ) ) );
}

void EvtSVVCP::decay( EvtParticle* p )
{
    //added by Lange Jan4,2000
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );

    double t;
    EvtId other_b;

    EvtCPUtil::getInstance()->OtherB( p, t, other_b, 0.5 );

    EvtComplex G0P, G1P, G1M;

    G1P = EvtComplex( getArg( 3 ) * cos( getArg( 4 ) ),
                      getArg( 3 ) * sin( getArg( 4 ) ) );
    G0P = EvtComplex( getArg( 5 ) * cos( getArg( 6 ) ),
                      getArg( 5 ) * sin( getArg( 6 ) ) );
    G1M = EvtComplex( getArg( 7 ) * cos( getArg( 8 ) ),
                      getArg( 7 ) * sin( getArg( 8 ) ) );

    EvtComplex lambda_km = EvtComplex( cos( -2 * getArg( 0 ) ),
                                       sin( -2 * getArg( 0 ) ) );

    double cdmt = cos( getArg( 1 ) * t / ( 2 * EvtConst::c ) );
    double sdmt = sin( getArg( 1 ) * t / ( 2 * EvtConst::c ) );

    EvtComplex cG0P, cG1P, cG1M;

    if ( other_b == B0B ) {
        cG0P = G0P * ( cdmt + lambda_km * EvtComplex( 0.0, getArg( 2 ) * sdmt ) );
        cG1P = G1P * ( cdmt + lambda_km * EvtComplex( 0.0, getArg( 2 ) * sdmt ) );
        cG1M = G1M * ( cdmt - lambda_km * EvtComplex( 0.0, getArg( 2 ) * sdmt ) );
    }
    if ( other_b == B0 ) {
        cG0P = G0P * ( cdmt + ( 1.0 / lambda_km ) *
                                  EvtComplex( 0.0, getArg( 2 ) * sdmt ) );
        cG1P = G1P * ( cdmt + ( 1.0 / lambda_km ) *
                                  EvtComplex( 0.0, getArg( 2 ) * sdmt ) );
        cG1M = -G1M * ( cdmt - ( 1.0 / lambda_km ) *
                                   EvtComplex( 0.0, getArg( 2 ) * sdmt ) );
    }

    EvtComplex A0, AP, AM;

    A0 = cG0P / sqrt( 2.0 );
    AP = ( cG1P + cG1M ) / sqrt( 2.0 );
    AM = ( cG1P - cG1M ) / sqrt( 2.0 );

    EvtSVVHelAmp::SVVHel( p, _amp2, getDaug( 0 ), getDaug( 1 ), AP, A0, AM );

    return;
}

std::string EvtSVVCP::getParamName( int i )
{
    switch ( i ) {
        case 0:
            return "weakPhase";
        case 1:
            return "deltaM";
        case 2:
            return "eta";
        case 3:
            return "G1Plus";
        case 4:
            return "G1PlusPhase";
        case 5:
            return "G0Plus";
        case 6:
            return "G0PlusPhase";
        case 7:
            return "G1Minus";
        case 8:
            return "G1MinusPhase";
        default:
            return "";
    }
}

std::string EvtSVVCP::getParamDefault( int i )
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
        case 7:
            return "1.0";
        case 8:
            return "0.0";
        default:
            return "";
    }
}
