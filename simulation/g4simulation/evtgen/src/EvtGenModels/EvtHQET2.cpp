
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

#include "EvtGenModels/EvtHQET2.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"

#include "EvtGenModels/EvtHQET2FF.hh"

#include <assert.h>
#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtHQET2::getName()
{
    return "HQET2";
}

EvtDecayBase* EvtHQET2::clone()
{
    return new EvtHQET2;
}

void EvtHQET2::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    calcamp->CalcAmp( p, _amp2, hqetffmodel.get() );
}

void EvtHQET2::initProbMax()
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

void EvtHQET2::init()
{
    checkNDaug( 3 );

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino
    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );
    if ( d1type == EvtSpinType::SCALAR ) {
        if ( getNArg() == 2 ) {
            hqetffmodel = std::make_unique<EvtHQET2FF>( getArg( 0 ), getArg( 1 ) );
            calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();

        } else if ( getNArg() == 3 ) {
            hqetffmodel = std::make_unique<EvtHQET2FF>( getArg( 0 ), getArg( 1 ),
                                                        getArg( 2 ) );
            calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();

        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "HQET2 model for scalar meson daughters needs 2 arguments for normal mode or 3 for extended. Sorry."
                << endl;
            ::abort();
        }

    } else if ( d1type == EvtSpinType::VECTOR ) {
        if ( getNArg() == 4 ) {
            hqetffmodel = std::make_unique<EvtHQET2FF>(
                getArg( 0 ), getArg( 1 ), getArg( 2 ), getArg( 3 ) );
            calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();

        } else if ( getNArg() == 5 ) {
            hqetffmodel = std::make_unique<EvtHQET2FF>( getArg( 0 ), getArg( 1 ),
                                                        getArg( 2 ), getArg( 3 ),
                                                        getArg( 4 ) );
            calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();

        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "HQET2 model for vector meson daughtersneeds 4 arguments for normal mode or 5 for extended. Sorry."
                << endl;
            ::abort();
        }

    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "HQET2 model handles only scalar and vector meson daughters. Sorry."
            << endl;
        ::abort();
    }
}
