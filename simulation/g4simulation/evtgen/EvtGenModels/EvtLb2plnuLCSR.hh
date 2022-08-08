
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

#ifndef EVTLB2PMUNULCSR_HH
#define EVTLB2PMUNULCSR_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

#include "EvtGenModels/EvtSLBaryonAmp.hh"

class EvtParticle;

// Description:Implementation of the Lb2plnuLCSR model
// Class to handle semileptonic Lb -> p l nu decays using the using form factor predictions from Light Cone Sum Rules.
// Description: Routine to implement Lb->p l nu semileptonic decays using form factor
//              predictions form Light Cone Sum Rules.  The form factors are from
//              A. Khodjamirian, C. Klein, T. Mannel and Y.-M. Wang, arXiv.1108.2971 (2011)

class EvtLb2plnuLCSR : public EvtDecayAmp {
  public:
    EvtLb2plnuLCSR();
    ~EvtLb2plnuLCSR();

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
