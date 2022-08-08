
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

#ifndef EVTBTOXSLL_HH
#define EVTBTOXSLL_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtBtoXsllUtil.hh"

#include <memory>

class EvtBtoXsllUtil;

// Description:
// Class to generate inclusive non-resonant B -> Xs l+ l- decays.
// Description: Routine to generate non-resonant B -> Xs l+ l- decays.
// It generates a dilepton mass spectrum according to Kruger and Sehgal
// and then generates the two lepton momenta accoring to Ali et al.
// The resultant X_s particles may be decayed by JETSET.

class EvtBtoXsll : public EvtDecayIncoherent {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* p ) override;

  private:
    std::unique_ptr<EvtBtoXsllUtil> _calcprob;
    double _dGdsProbMax;
    double _dGdsdupProbMax;
    double _mb;
    double _ms;
    double _mq;
    double _pf;
    double _mxmin;
};

#endif
