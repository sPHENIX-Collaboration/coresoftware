
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

#ifndef EVTRARELBTOLLL_HH
#define EVTRARELBTOLLL_HH 1

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"
#include "EvtGenModels/EvtRareLbToLllWC.hh"

#include <memory>

// Description:
//      Implements the rare Lb --> Lambda^(*) ell ell models described in
//      http://arxiv.org/pdf/1108.6129.pdf

class EvtRareLbToLll : public EvtDecayAmp {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void init() override;

    void initProbMax() override;

    void decay( EvtParticle* parent ) override;

  protected:
    void calcAmp( EvtAmp& amp, EvtParticle* parent );

    void HadronicAmp( EvtParticle* parent, EvtParticle* lambda, EvtVector4C* T,
                      const int i, const int j );

    void HadronicAmpRS( EvtParticle* parent, EvtParticle* lambda,
                        EvtVector4C* T, const int i, const int j );

    bool isParticle( EvtParticle* parent ) const;

  private:
    double m_maxProbability;

    std::unique_ptr<EvtRareLbToLllFFBase> ffmodel_;
    std::unique_ptr<EvtRareLbToLllWC> wcmodel_;
};
#endif    //
