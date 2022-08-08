
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

#include "EvtGenModels/EvtHQET.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"

#include "EvtGenModels/EvtHQETFF.hh"

#include <assert.h>
#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtHQET::getName()
{
    return "HQET";
}

EvtDecayBase* EvtHQET::clone()
{
    return new EvtHQET;
}

void EvtHQET::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    calcamp->CalcAmp( p, _amp2, hqetffmodel.get() );
}

void EvtHQET::initProbMax()
{
    EvtId parnum, mesnum, lnum, nunum;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    lnum = getDaug( 1 );
    nunum = getDaug( 2 );

    double mymaxprob = calcamp->CalcMaxProb( parnum, mesnum, lnum, nunum,
                                             hqetffmodel.get() );

    setProbMax( mymaxprob );
}

void EvtHQET::init()
{
    checkNDaug( 3 );

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino
    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );
    if ( d1type == EvtSpinType::SCALAR ) {
        checkNArg( 1, 2 );
        if ( getNArg() == 1 )
            hqetffmodel = std::make_unique<EvtHQETFF>( getArg( 0 ) );
        else
            hqetffmodel = std::make_unique<EvtHQETFF>( getArg( 0 ), getArg( 1 ) );
        calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();
    } else if ( d1type == EvtSpinType::VECTOR ) {
        checkNArg( 3, 4 );
        if ( getNArg() == 3 )
            hqetffmodel = std::make_unique<EvtHQETFF>( getArg( 0 ), getArg( 1 ),
                                                       getArg( 2 ) );
        else
            hqetffmodel = std::make_unique<EvtHQETFF>( getArg( 0 ), getArg( 1 ),
                                                       getArg( 2 ), getArg( 3 ) );
        calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "HQET model handles only scalar and vector meson daughters. Sorry."
            << endl;
        ::abort();
    }
}
