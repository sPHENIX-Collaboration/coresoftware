
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

#ifndef EVTPROPSLPOLE_HH
#define EVTPROPSLPOLE_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtPoint1D.hh"
#include "EvtGenBase/EvtSemiLeptonicAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

#include <memory>

class Evtparticle;

// Description:Semileptonic decays with pole form form factors

class EvtPropSLPole : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void initProbMax() override;
    void init() override;

    double calBreitWigner( EvtParticle* pmeson, EvtPoint1D point );
    double calBreitWignerBasic( double maxMass );

    double calcMaxProb( EvtId parent, EvtId meson, EvtId lepton, EvtId nudaug,
                        EvtSemiLeptonicFF* FormFactors );

  private:
    bool _includeDecayFact;
    bool _includeBirthFact;
    double _mass;
    double _massMin;
    double _massMax;
    double _width;
    double _maxRange;
    EvtSpinType::spintype _spin;

    double _blatt;
    bool _isProbMaxSet;

    std::unique_ptr<EvtSemiLeptonicFF> SLPoleffmodel;
    std::unique_ptr<EvtSemiLeptonicAmp> calcamp;
};

#endif
