
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

#ifndef EVTNORADCORR_HH
#define EVTNORADCORR_HH

#include "EvtGenBase/EvtAbsRadCorr.hh"

#include <string>

class EvtParticle;

// Description: Create an empty radiative correction engine which does nothing.
// This is required since the EvtGen constructor still needs at least one
// concrete implementation of EvtAbsRadCorr for the case when Photos is not used.

class EvtNoRadCorr : public EvtAbsRadCorr {
  public:
    EvtNoRadCorr() { ; }
    virtual ~EvtNoRadCorr() { ; }

    virtual void doRadCorr( EvtParticle* ) { ; }

  private:
};

#endif
