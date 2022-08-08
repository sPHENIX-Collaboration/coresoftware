
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

#ifndef EVTABS_EXTERNALGEN_HH
#define EVTABS_EXTERNALGEN_HH

#include "EvtGenBase/EvtParticle.hh"

// Description: Pure abstract interface for external physics generators

class EvtAbsExternalGen {
  public:
    virtual ~EvtAbsExternalGen() = default;
    virtual bool doDecay( EvtParticle* theMother ) = 0;
    virtual double getDecayProb( EvtParticle* ) { return 1.0; }
    virtual void initialise() = 0;
};

#endif
