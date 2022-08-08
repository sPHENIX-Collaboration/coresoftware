
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

#ifndef EVTBSTOLLGAMMAISRFSR_AMP_HH
#define EVTBSTOLLGAMMAISRFSR_AMP_HH

class EvtId;
class EvtAmp;
class EvtParticle;
class Evtbs2llGammaFF;
class EvtbTosllWilsCoeffNLO;

// Description: Preparation of the decay amplitude for the process:
//              B^0_{q}(p,M1) -> gamma(k) ell^+(p1,m) ell^-(p2,m).
//              See the Internal LHCb Note LHCb-INT-2011-011.
//
// Note: The code of this module is based on the EvtbTosllVectorAmp.cpp
//	 module code.
//	 The main functiom for the amplitude calculation retuns the
//	 amplitude for the decay  B -> gamma ell^+ ell^-
//	 In our calculations we assume, that photon is the first
//       daughter particle (iG=0) and leptons are the second and thirds
//       daughter particles (il1=1 and il2=2).

class Evtbs2llGammaISRFSRAmp {
  public:
    void CalcAmp( EvtParticle* parent, EvtAmp& amp, Evtbs2llGammaFF* formFactors,
                  EvtbTosllWilsCoeffNLO* WilsCoeff, double mu, int Nf, int sr,
                  int res_swch, int ias, double Egamma_min, double CKM_A,
                  double CKM_lambda, double CKM_barrho, double CKM_bareta,
                  double mumumass_min );

    double CalcMaxProb( EvtId parnum, EvtId photnum, EvtId l1num, EvtId l2num,
                        Evtbs2llGammaFF* formFactors,
                        EvtbTosllWilsCoeffNLO* WilsCoeff, double mu, int Nf,
                        int sr, int res_swch, int ias, double Egamma_min,
                        double CKM_A, double CKM_lambda, double CKM_barrho,
                        double CKM_bareta, double mumumass_min );

    double lambda( double a, double b, double c );
};

#endif
