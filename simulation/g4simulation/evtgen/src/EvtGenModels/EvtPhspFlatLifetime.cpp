
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

#include "EvtGenModels/EvtPhspFlatLifetime.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

//==============================================================================
// Return the name of the model
//==============================================================================
std::string EvtPhspFlatLifetime::getName()
{
    return "PHSPFLATLIFETIME";
}

//==============================================================================
// Copy the model
//==============================================================================
EvtDecayBase* EvtPhspFlatLifetime::clone()
{
    return new EvtPhspFlatLifetime;
}

//==============================================================================
// Initialize the model
//==============================================================================
void EvtPhspFlatLifetime::init()
{
    // check that there is 1 argument in the decay file
    checkNArg( 1 );
    // this argument is the lifetime upper edge (in ps)
    m_maxLifetime = getArg( 0 ) * EvtConst::c * 1.e-12;
}

//==============================================================================
// Compute the maximum probability (max of the pdf)
//==============================================================================
void EvtPhspFlatLifetime::initProbMax()
{
    noProbMax();
}

//==============================================================================
// Decay the particle according to the model
//==============================================================================
void EvtPhspFlatLifetime::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    // generate the lifetime flat between 0 and max
    double l = EvtRandom::Flat( 0., m_maxLifetime );
    // modify the lifetime of the particle (in mm)
    p->setLifetime( l );
}
