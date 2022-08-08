
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

#ifndef EVTBSTOLLLL_AMP_HH
#define EVTBSTOLLLL_AMP_HH

class EvtId;
class EvtAmp;
class EvtParticle;
class Evtbs2llGammaFF;
class EvtbTosllWilsCoeffNLO;

class EvtbsToLLLLAmp {
  public:
    void CalcAmp( EvtParticle* parent, EvtAmp& amp, Evtbs2llGammaFF* formFactors,
                  EvtbTosllWilsCoeffNLO* WilsCoeff, double mu, int Nf,
                  int res_swch, int ias, double CKM_A, double CKM_lambda,
                  double CKM_barrho, double CKM_bareta );

    double CalcMaxProb(
        //                       EvtId parnum,
        //                       EvtId l1num, EvtId l2num,
        //                       EvtId l3num, EvtId l4num,
        //		         Evtbs2llGammaFF *formFactors,
        //                       EvtbTosllWilsCoeffNLO *WilsCoeff,
        //                       double mu, int Nf, int res_swch, int ias,
        //                       double CKM_A, double CKM_lambda,
        //                       double CKM_barrho, double CKM_bareta
    );

    double lambda( double a, double b, double c );
};

#endif
