
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

#ifndef EVTBLLNUL_HH
#define EVTBLLNUL_HH

#include "EvtGenBase/EvtDecayAmp.hh"

#include "EvtGenModels/EvtBLLNuLAmp.hh"

#include <string>

class EvtParticle;
class EvtbTosllMSFF;    // Form factor class

// Description: The header file for the model "BLLNUL" which simulates
//              the rare four-leptonic B-decays
//              B^-(p) -> ell^+(k_1) ell^-(k_2) neu (k_3) ell^-(k_4)

class EvtBLLNuL : public EvtDecayAmp {
  public:
    EvtBLLNuL();

    virtual std::string getName() override;
    virtual EvtDecayBase* clone() override;

    virtual void init() override;
    virtual void initProbMax() override;
    virtual void decay( EvtParticle* p ) override;

  private:
    EvtBLLNuLAmp calcAmp_;
};

#endif
