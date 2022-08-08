
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

#include "EvtGenModels/EvtBcPsiNPi.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtBcPsiNPi::EvtBcPsiNPi()
{
    _beta = -0.108;
    _mRho = 0.775;
    _gammaRho = 0.149;
    _mRhopr = 1.364;
    _gammaRhopr = 0.400;
    _mA1 = 1.23;
    _gammaA1 = 0.4;

    FA0_N = 5.9;
    FA0_c1 = 0.049;
    FA0_c2 = 0.0015;
    FAm_N = 0.0;
    FAm_c1 = 0.0;
    FAm_c2 = 0.0;
    FAp_N = -0.074;
    FAp_c1 = 0.049;
    FAp_c2 = 0.0015;
    FV_N = 0.11;
    FV_c1 = 0.049;
    FV_c2 = 0.0015;
}

std::string EvtBcPsiNPi::getName()
{
    return "BC_PSI_NPI";
}

EvtBcPsiNPi* EvtBcPsiNPi::clone()
{
    return new EvtBcPsiNPi;
}

void EvtBcPsiNPi::init()
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

void EvtBcPsiNPi::initProbMax()
{
    setProbMax( 100. );
    if ( getNDaug() == 2 ) {
        setProbMax( 330. );
    } else if ( getNDaug() == 3 ) {
        setProbMax( 11000. );    // checked with 30k events
    } else if ( getNDaug() == 4 ) {
        setProbMax( 36000. );
    }
}
