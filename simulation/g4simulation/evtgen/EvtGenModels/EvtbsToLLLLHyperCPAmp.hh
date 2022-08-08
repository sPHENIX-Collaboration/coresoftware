
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

#ifndef EVTBSTOLLLL_HYPERCPAMP_HH
#define EVTBSTOLLLL_HYPERCPAMP_HH

class EvtId;
class EvtAmp;
class EvtParticle;

// Description: Preparation of the decay amplitude for the process:
//              B^0_{q}(p,M1) -> ell^+(k1,m) ell^-(k2,m) ell^+(k3,m) ell^-(k4,m)
//              in the HyperCP model.
//
//              [1] D.S.Gorbunov, Nucl.Phys.B602, pp.213-237 (2001);
//              [2] S.V. Demidov, D.S.Gorbunov, hep-ph/1112.5230v2, 17 April 2012.
//
// Note: The code of this module is based on the EvtbsToLLLLAmp.cpp module code.

class EvtbsToLLLLHyperCPAmp {
  public:
    void CalcAmp( EvtParticle* parent, EvtAmp& amp, double mS, double mP,
                  double gammaS, double gammaP, double mLiiLR, double Fc,
                  double mD23LL, double mD23RR, double mD32LL, double mD32RR,
                  double mD13LL, double mD13RR, double mD31LL, double mD31RR );

    double CalcMaxProb( EvtId parnum, EvtId l1num, EvtId l2num, EvtId l3num,
                        EvtId l4num, double mS, double mP, double gammaS,
                        double gammaP, double mLiiLR, double Fc, double mD23LL,
                        double mD23RR, double mD32LL, double mD32RR,
                        double mD13LL, double mD13RR, double mD31LL,
                        double mD31RR );

    double lambda( double a, double b, double c );
};

#endif
