
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

#include "EvtGenModels/EvtPhsp.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

std::string EvtPhsp::getName()
{
    return "PHSP";
}

EvtDecayBase* EvtPhsp::clone()
{
    return new EvtPhsp;
}

void EvtPhsp::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );
}

void EvtPhsp::initProbMax()
{
    noProbMax();
}

void EvtPhsp::decay( EvtParticle* p )
{
    //unneeded - lange - may13-02
    //if ( p->getNDaug() != 0 ) {
    //Will end up here because maxrate multiplies by 1.2
    //  EvtGenReport(EVTGEN_DEBUG,"EvtGen") << "In EvtPhsp: has "
    //			   <<" daugthers should not be here!"<<endl;
    //  return;
    //}

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    return;
}
