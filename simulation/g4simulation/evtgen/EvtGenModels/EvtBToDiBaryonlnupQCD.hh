
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

#ifndef EVTBTODIBARYONLNUPQCD_HH
#define EVTBTODIBARYONLNUPQCD_HH

#include "EvtGenBase/EvtDecayAmp.hh"

#include "EvtGenModels/EvtBToDiBaryonlnupQCDFF.hh"
#include "EvtGenModels/EvtSLDiBaryonAmp.hh"

#include <memory>
#include <string>

class EvtParticle;

// Description: Class to handle semileptonic B -> Baryon Anti-baryon l nu decays
// using the using form factor predictions from pQCD counting rules. Taken
// from arXiv:1107.0801

class EvtBToDiBaryonlnupQCD : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void initProbMax() override;
    void init() override;

  private:
    std::unique_ptr<EvtBToDiBaryonlnupQCDFF> ffModel_;
    std::unique_ptr<EvtSLDiBaryonAmp> calcAmp_;
};

#endif
