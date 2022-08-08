
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

#include "EvtGenModels/EvtBToVlnuBall.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"

#include "EvtGenModels/EvtBToVlnuBallFF.hh"

#include <assert.h>
#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtBToVlnuBall::getName()
{
    return "BTOVLNUBALL";
}

EvtBToVlnuBall* EvtBToVlnuBall::clone()
{
    return new EvtBToVlnuBall;
}

void EvtBToVlnuBall::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    _calcamp->CalcAmp( p, _amp2, _Ballmodel.get() );
}

void EvtBToVlnuBall::initProbMax()
{
    EvtId parnum, mesnum, lnum, nunum;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    lnum = getDaug( 1 );
    nunum = getDaug( 2 );

    double mymaxprob = _calcamp->CalcMaxProb( parnum, mesnum, lnum, nunum,
                                              _Ballmodel.get() );

    setProbMax( mymaxprob );
}

void EvtBToVlnuBall::init()
{
    checkNDaug( 3 );

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino
    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );
    if ( d1type == EvtSpinType::VECTOR ) {
        checkNArg( 8 );    // the number of arguments needed for the Ball model
        _Ballmodel = std::make_unique<EvtBToVlnuBallFF>(
            getArg( 0 ), getArg( 1 ), getArg( 2 ), getArg( 3 ), getArg( 4 ),
            getArg( 5 ), getArg( 6 ), getArg( 7 ) );
        _calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Ball model handles only vector meson daughters. Sorry." << endl;
        ::abort();
    }
}
