
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

#ifndef EVTBS2LLGAMMAMNT_HH
#define EVTBS2LLGAMMAMNT_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;
class Evtbs2llGammaFF;    // my class with ff for rare semileptonic B-decays
class Evtbs2llGammaAmp;    // my class with amplitudes for rare semileptonic B-decays
class EvtbTosllWilsCoeffNLO;    // my class with Wilson coefficients in NLO

// Description: the B -> Gamma ell^+ ell^- decay channel description
//		based on the papers:
//		D.Melikhov, N.Nikitin, K.Toms, Yad. Fiz. 62, No 11
//              and
//              I.Balakireva, D.Melikhov, N.Nikitin, D.Tlisov,
//                                           e-Print: arXiv:0911.0605 [hep-ph].

class Evtbs2llGammaMNT : public EvtDecayAmp {
  public:
    Evtbs2llGammaMNT() {}
    virtual ~Evtbs2llGammaMNT();

    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  private:
    Evtbs2llGammaFF* _mntffmodel;
    Evtbs2llGammaAmp* _calcamp;
    EvtbTosllWilsCoeffNLO* _wilscoeff;
};

#endif
