
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

#ifndef EVTBTOSLLMS_HH
#define EVTBTOSLLMS_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;
class EvtbTosllFFNew;    // my class with ff for rare semileptonic B-decays
class EvtbTosllAmpNew;    // my class with amplitudes for rare semileptonic B-decays
class EvtbTosllWilsCoeffNLO;    // my class with Wilson coefficients in NLO

// Description: Implementation of the B -> K* ell^+ ell^- decay according
//              to the paper: D.Melikhov, B.Stech, PRD62, 014006 (2000).
//
//              This is the derived class of the base class "EvtDecayAmp",
//              but the decay amplitude will be calculated in the class
//              "EvtbTosllAmpNew" (see file EvtbTosllAmpNew.hh).

class EvtbTosllMS : public EvtDecayAmp {
  public:
    EvtbTosllMS(){};
    virtual ~EvtbTosllMS();

    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  private:
    EvtbTosllFFNew* _msffmodel;
    EvtbTosllAmpNew* _calcamp;
    EvtbTosllWilsCoeffNLO* _wilscoeff;
};

#endif
