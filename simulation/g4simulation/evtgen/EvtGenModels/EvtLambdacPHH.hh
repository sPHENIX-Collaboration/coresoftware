
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

#ifndef EVTLAMBDACPHH_HH
#define EVTLAMBDACPHH_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <string>
#include <vector>

class EvtParticle;

// Description: Decay model for Lambda_c -> K- pi+ p using amplitudes
//              from the Fermilab E791 analysis: arXiv:hep-ex/9912003v1

class EvtLambdacPHH : public EvtDecayAmp {
  public:
    EvtLambdacPHH();

    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  protected:
    // Resonance enumeration
    enum LcResLabel
    {
        NonReson = 0,
        Kstar,
        Delta,
        Lambda
    };

    // Amplitude functions
    std::vector<EvtComplex> calcResAmpTerms( EvtLambdacPHH::LcResLabel resIndex,
                                             const EvtResonance2& res,
                                             double norm ) const;

    EvtComplex DecayAmp3( EvtLambdacPHH::LcResLabel resonance, int m,
                          int mprime, double theta_res, double phi_res,
                          double theta_prime_daughter_res,
                          double phi_prime_daughter_res ) const;

    EvtComplex fampl3( double amplitude_res, double phi_res, int spinMother,
                       int m_spinMother, int m_prime_spinMother,
                       double theta_res, float spin_res, float m_spin_res,
                       float m_prime_spin_res, double theta_daughter_res,
                       double phi_prime_daughter_res ) const;

    // Find resonance normalisation factors
    void calcNormalisations();

    void getFitFractions();

    // Inverse cos/sin functions that checks for valid arguments
    double getACos( double num, double denom ) const;
    double getASin( double num, double denom ) const;

  private:
    // Daughter ordering for K-, pi+, p
    int _d1, _d2, _d3;

    // Resonance parameters
    double _Nplusplus, _Nplusminus, _Nminusplus, _Nminusminus;
    double _phiNplusplus, _phiNplusminus, _phiNminusplus, _phiNminusminus;
    double _E1, _phiE1, _E2, _phiE2, _E3, _phiE3, _E4, _phiE4;
    double _F1, _phiF1, _F2, _phiF2, _H1, _phiH1, _H2, _phiH2;

    double _NRNorm, _KstarNorm, _DeltaNorm, _LambdaNorm;
    double _KstarM, _KstarW, _KstarR;
    double _DeltaM, _DeltaW, _DeltaR;
    double _LambdaM, _LambdaW, _LambdaR;
    double _Lambda_cR;

    EvtVector4R _zprime, _p4_Lambda_c;
    double _zpMag, _p4_Lambdac_Mag;
};

#endif
