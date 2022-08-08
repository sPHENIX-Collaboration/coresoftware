
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

#include "EvtGenModels/EvtBcBsNPi.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtBcBsNPi::EvtBcBsNPi()
{
    _beta = -0.108;
    _mRho = 0.775;
    _gammaRho = 0.149;
    _mRhopr = 1.364;
    _gammaRhopr = 0.400;
    _mA1 = 1.23;
    _gammaA1 = 0.4;
    //		Fp_N=1.3; Fp_c1=0.30; Fp_c2=0.069;
    Fp_N = 3 * 1.3;
    Fp_c1 = 0.30;
    Fp_c2 = 0.069;
    Fm_N = 0.0;
    Fm_c1 = 0.0;
    Fm_c2 = 0.0;
}

std::string EvtBcBsNPi::getName()
{
    return "BC_BS_NPI";
}

EvtBcBsNPi* EvtBcBsNPi::clone()
{
    return new EvtBcBsNPi;
}

void EvtBcBsNPi::init()
{
    checkNArg( 0 );

    // check spins
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    // the others are scalar
    for ( int i = 1; i <= ( getNDaug() - 1 ); i++ ) {
        checkSpinDaughter( i, EvtSpinType::SCALAR );
    }
}

void EvtBcBsNPi::initProbMax()
{
    if ( getNDaug() == 2 ) {
        setProbMax( 250. );
    } else if ( getNDaug() == 3 ) {
        setProbMax( 25000. );    // checked at 30k events
    } else if ( getNDaug() == 4 ) {
        setProbMax( 45000. );    // checked at 30k events
    }
}
