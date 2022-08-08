
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

#ifndef EvtBcToNPi_HH
#define EvtBcToNPi_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <string>

// Description: General decay model for Bc -> V + npi and Bc -> P + npi

class EvtBcToNPi : public EvtDecayAmp {
  public:
    EvtBcToNPi( bool printAuthorInfo = false );

    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* p ) override;

  protected:
    int nCall;
    double maxAmp2;

    // Bc form factors
    double _maxProb;
    double FA0_N, FA0_c1, FA0_c2;
    double FAm_N, FAm_c1, FAm_c2;
    double FAp_N, FAp_c1, FAp_c2;
    double FV_N, FV_c1, FV_c2;

    double Fp_N, Fp_c1, Fp_c2;
    double Fm_N, Fm_c1, Fm_c2;

    // W -> pi... form factors
    double _beta;
    double _mRho;
    double _gammaRho;
    double _mRhopr;
    double _gammaRhopr;
    double _mA1;
    double _gammaA1;

    double _ee( double M, double m1, double m2 );
    double _pp( double M, double m1, double m2 );
    EvtComplex Fpi( EvtVector4R q1, EvtVector4R q2 );
    double pi3G( double m2, int dupD );

  private:
    void printAuthorInfo();
};

#endif
