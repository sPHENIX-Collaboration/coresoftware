
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

#include "EvtGenModels/EvtBcBsStarNPi.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtBcBsStarNPi::EvtBcBsStarNPi()
{
    _beta = -0.108;
    _mRho = 0.775;
    _gammaRho = 0.149;
    _mRhopr = 1.364;
    _gammaRhopr = 0.400;
    _mA1 = 1.23;
    _gammaA1 = 0.4;

    FA0_N = 8.1;
    FA0_c1 = 0.30;
    FA0_c2 = 0.069;
    FAm_N = 0.0;
    FAm_c1 = 0.0;
    FAm_c2 = 0.0;
    FAp_N = 0.15;
    FAp_c1 = 0.30;
    FAp_c2 = 0.069;
    FV_N = 1.08;
    FV_c1 = 0.30;
    FV_c2 = 0.069;
}

std::string EvtBcBsStarNPi::getName()
{
    return "BC_BSSTAR_NPI";
}

EvtBcBsStarNPi* EvtBcBsStarNPi::clone()
{
    return new EvtBcBsStarNPi;
}

void EvtBcBsStarNPi::init()
{
    checkNArg( 0 );

    // check spins
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    // the others are scalar
    for ( int i = 1; i <= ( getNDaug() - 1 ); i++ ) {
        checkSpinDaughter( i, EvtSpinType::SCALAR );
    }
}

void EvtBcBsStarNPi::initProbMax()
{
    if ( getNDaug() == 2 ) {
        setProbMax( 100. );
    } else if ( getNDaug() == 3 ) {
        setProbMax( 40000. );
    } else if ( getNDaug() == 4 ) {
        setProbMax( 620. );    // checked, 30k events
    }
}
