
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

#include "EvtGenExternal/EvtTauola.hh"

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenModels/EvtAbsExternalGen.hh"

#include "EvtGenExternal/EvtExternalGenFactory.hh"

#include <cmath>
#include <iostream>
#include <string>

std::string EvtTauola::getName()
{
    return "TAUOLA";
}

EvtDecayBase* EvtTauola::clone()
{
    return new EvtTauola();
}

void EvtTauola::init()
{
}

void EvtTauola::initProbMax()
{
    noProbMax();
}

void EvtTauola::decay( EvtParticle* p )
{
    // We check to see if the Tauola engine has been created before doing the decay.
    // This should only create the full Tauola engine once, and all clones will
    // point to the same engine.

    if ( !_tauolaEngine ) {
        _tauolaEngine = EvtExternalGenFactory::getInstance()->getGenerator(
            EvtExternalGenFactory::TauolaGenId );
    }

    if ( _tauolaEngine ) {
        _tauolaEngine->doDecay( p );
    }
}
