
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

#ifndef EvtBcVHad_HH
#define EvtBcVHad_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtBCVFF2.hh"
#include "EvtGenModels/EvtWHad.hh"

#include <memory>
#include <string>

class EvtParticle;

// Description: Module to implement Bc -> psi + (n pi) + (m K) decays

class EvtBcVHad : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;
    void initProbMax() override;
    void init() override;
    void decay( EvtParticle* p ) override;

  protected:
    // Hadronic current function
    EvtVector4C hardCurr( EvtParticle* root_particle ) const;

  private:
    // whichfit --- code of the Bc -> VW formfactor set:
    //   1 - SR
    //   2 - PM
    int whichfit;

    // idVector --- final vector particle code
    int idVector;

    // out_code: code of the hadronic final state
    //   1 - pi+
    //   2 - pi+ pi0
    //   3 -  pi+ pi+ pi-
    //   4 - 4pi
    //   5 - pi+ pi+ pi- pi- pi+
    //   6 - K+ K- pi+
    //   7 - K+ pi+ pi-
    //   8 - K_S0 K+
    int out_code;

    std::unique_ptr<EvtBCVFF2> ffmodel;
    std::unique_ptr<EvtWHad> wcurr;
};

#endif
