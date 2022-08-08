
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

#include "EvtGenModels/EvtKKLambdaC.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicBaryonAmp.hh"

#include "EvtGenModels/EvtKKLambdaCFF.hh"

#include <stdlib.h>
#include <string>

std::string EvtKKLambdaC::getName()
{
    return "KK_LAMBDAC_SL";
}

EvtDecayBase* EvtKKLambdaC::clone()
{
    return new EvtKKLambdaC;
}

void EvtKKLambdaC::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    _calcamp->CalcAmp( p, _amp2, _ffmodel.get() );
}

void EvtKKLambdaC::initProbMax()
{
    EvtId parnum, mesnum, lnum, nunum;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    lnum = getDaug( 1 );
    nunum = getDaug( 2 );

    //double mymaxprob = _calcamp->CalcMaxProb(parnum,mesnum,
    //                           lnum,nunum,_ffmodel);
    double mymaxprob = 100.;
    setProbMax( mymaxprob );
}

void EvtKKLambdaC::init()
{
    checkNDaug( 3 );

    //We expect the parent to be a dirac
    //and the daughters to be dirac lepton neutrino

    checkSpinParent( EvtSpinType::DIRAC );
    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    _ffmodel = std::make_unique<EvtKKLambdaCFF>( getNArg(), getArgs() );

    _calcamp = std::make_unique<EvtSemiLeptonicBaryonAmp>();
}
