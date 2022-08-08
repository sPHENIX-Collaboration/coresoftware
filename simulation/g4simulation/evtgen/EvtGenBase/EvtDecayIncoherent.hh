
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

#ifndef EvtDecayIncoherent_HH
#define EvtDecayIncoherent_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"

// Description: Base class for models that calculate
//              decay kinematics and do not do any accept/reject.
//              Useful e.g. for interface to other generators

class EvtDecayIncoherent : public EvtDecayBase {
  public:
    void makeDecay( EvtParticle* p, bool recursive = true ) override;

    virtual ~EvtDecayIncoherent() {}

    void setDaughterSpinDensity( int daughter )
    {
        spinDensitySet[daughter] = 1;
        return;
    }

    int isDaughterSpinDensitySet( int daughter )
    {
        return spinDensitySet[daughter];
    }

  private:
    int spinDensitySet[MAX_DAUG];
};

#endif
