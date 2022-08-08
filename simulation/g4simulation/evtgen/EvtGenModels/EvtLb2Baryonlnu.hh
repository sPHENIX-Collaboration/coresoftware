
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

#ifndef EVTLB2BARYONLNU_HH
#define EVTLB2BARYONLNU_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

#include "EvtGenModels/EvtSLBaryonAmp.hh"

class EvtParticle;

// Description:Implementation of the Lb2Baryonlnu model
// Class to handle semileptonic Lb -> p l nu decays using the using form factor predictions based on the quark model.  Here baryon can be Lc+, p, Lc*+, N*+.
// Description: Routine to implement Lb->N*+ l nu semileptonic decays using form factor predictions based on the quark model.  The form factors are from W. Roberts, M. Pervin, S. Chapstick, (2011). arXiv:nucl-th/0503030v1.  The model can be used for decays to all N*+ states with J^{P} = 1/2^{+}, 1/2^{-}, 3/2^{+}, 3/2^{-} and, in addition, decays to p, Lc+, Lc(2593)+, Lc(2625)+

class EvtLb2Baryonlnu : public EvtDecayAmp {
  public:
    EvtLb2Baryonlnu();
    ~EvtLb2Baryonlnu();

    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void initProbMax() override;
    void init() override;

  private:
    EvtSemiLeptonicFF* ffmodel;
    EvtSLBaryonAmp* calcamp;
};

#endif
