
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

#include "EvtGenModels/EvtBToPlnuBK.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"

#include "EvtGenModels/EvtBToPlnuBKFF.hh"

#include <assert.h>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::fstream;

std::string EvtBToPlnuBK::getName()
{
    return "BTOPLNUBK";
}

EvtBToPlnuBK* EvtBToPlnuBK::clone()
{
    return new EvtBToPlnuBK;
}

void EvtBToPlnuBK::initProbMax()
{
    EvtId parnum, mesnum, lnum, nunum;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    lnum = getDaug( 1 );
    nunum = getDaug( 2 );

    double mymaxprob = calcamp->CalcMaxProb( parnum, mesnum, lnum, nunum,
                                             BKmodel.get() );

    setProbMax( mymaxprob );
}

void EvtBToPlnuBK::init()
{
    checkNDaug( 3 );

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino
    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );
    if ( d1type == EvtSpinType::SCALAR ) {
        checkNArg( 2 );
        BKmodel = std::make_unique<EvtBToPlnuBKFF>( getArg( 0 ), getArg( 1 ) );
        calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "BK model handles only scalar meson daughters. Sorry." << endl;
        ::abort();
    }
}

void EvtBToPlnuBK::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    calcamp->CalcAmp( p, _amp2, BKmodel.get() );
}
