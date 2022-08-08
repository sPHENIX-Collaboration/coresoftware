
/***********************************************************************
* Copyright 1998-2021 CERN for the benefit of the EvtGen authors       *
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

#ifndef EvtBTOXELNU_HH
#define EvtBTOXELNU_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

#include <memory>

class EvtParticle;

/** The class provides the form factors for orbitally excited semileptonic decays
 */
class EvtBToXElNu : public EvtDecayAmp {
  public:
    /** Default constructor */
    EvtBToXElNu() = default;

    /** Returns name of module */
    std::string getName() override;

    /** Clones module */
    EvtDecayBase* clone() override;

    /** Creates a decay */
    void decay( EvtParticle* p ) override;

    /** Sets maximal probab. */
    void initProbMax() override;

    /** Initializes module */
    void init() override;

  private:
    /** Pointers needed for FFs */
    std::unique_ptr<EvtSemiLeptonicFF> m_ffmodel{ nullptr };

    /** Pointers needed to calculate amplitude */
    std::unique_ptr<EvtSemiLeptonicAmp> m_calcamp{ nullptr };
};
#endif
