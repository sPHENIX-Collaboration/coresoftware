
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

#ifndef EVTPHSPFLATLIFETIME_HH
#define EVTPHSPFLATLIFETIME_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

class EvtParticle;

// Description:
//   Class to handle generic phase space decays not done
//   in other decay models, with a flat lifetime

class EvtPhspFlatLifetime : public EvtDecayIncoherent {
  public:
    /// Constructor
    EvtPhspFlatLifetime() : m_maxLifetime( 0. ){};

    /// Destructor
    virtual ~EvtPhspFlatLifetime(){};

    /// return name of the model
    std::string getName() override;

    /// Clone
    EvtDecayBase* clone() override;

    /// Compute maximum weight
    void initProbMax() override;

    /// Initialize the model
    void init() override;

    /// Perform the decay
    void decay( EvtParticle* p ) override;

  private:
    /// parameter of the model: maximum of the generated lifetime (in ps)
    double m_maxLifetime;
};

#endif
