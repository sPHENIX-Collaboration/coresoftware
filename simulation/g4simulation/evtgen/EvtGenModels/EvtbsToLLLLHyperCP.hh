
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

#ifndef EVTBSTOLLLL_HYPERCP_HH
#define EVTBSTOLLLL_HYPERCP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;
class EvtbsToLLLLHyperCPAmp;    // my class with amplitudes for rare four-leptonic B-decays

class EvtbsToLLLLHyperCP : public EvtDecayAmp {
  public:
    EvtbsToLLLLHyperCP(){};
    virtual ~EvtbsToLLLLHyperCP();

    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  private:
    EvtbsToLLLLHyperCPAmp* _calcamp;
};

#endif
