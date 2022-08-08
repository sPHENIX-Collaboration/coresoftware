
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

#include "EvtGenModels/EvtPhspDecaytimeCut.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>
#include <iostream>

std::string EvtPhspDecaytimeCut::getName()
{
    return "PHSPDECAYTIMECUT";
}

EvtDecayBase* EvtPhspDecaytimeCut::clone()
{
    return new EvtPhspDecaytimeCut;
}

void EvtPhspDecaytimeCut::init()
{
    // check that there are 1 arguments
    checkNArg( 1 );
    // This argument is minimum decay time in ps converted here to EvtGen
    // units in which c=1
    m_minDecayTime = getArg( 0 ) * EvtConst::c * 1.e-12;
}

void EvtPhspDecaytimeCut::initProbMax()
{
    noProbMax();
}

void EvtPhspDecaytimeCut::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    // Shift generated decay time by minimum we require
    const double currentDecaytime = p->getLifetime();
    p->setLifetime( currentDecaytime + m_minDecayTime );

    return;
}
