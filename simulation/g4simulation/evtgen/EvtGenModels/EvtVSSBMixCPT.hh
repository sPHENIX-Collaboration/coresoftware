
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

#ifndef EVTVSSBMIXCPT_HH
#define EVTVSSBMIXCPT_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtParticle.hh"

// Description:
//    Routine to decay vector-> scalar scalar with coherent BB-like mixing
//                              including CPT effects
//    Based on VSSBMIX

class EvtVSSBMixCPT : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void init() override;
    void initProbMax() override;

    int nRealDaughters() override { return 2; }

    std::string getParamName( int i ) override;
    std::string getParamDefault( int i ) override;

  private:
    double _freq;    // mixing frequency in hbar/mm
    double _dGamma;
    EvtComplex _qoverp;
    EvtComplex _poverq;
    EvtComplex _z;
    double _chib0_b0bar;
    double _chib0bar_b0;

    EvtComplex _A_f;
    EvtComplex _Abar_f;

    EvtComplex _A_fbar;
    EvtComplex _Abar_fbar;
};

#endif
