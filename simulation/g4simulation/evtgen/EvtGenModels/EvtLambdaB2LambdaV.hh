
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

#ifndef EVTLAMBDAB2LAMBDAV_HH
#define EVTLAMBDAB2LAMBDAV_HH

#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

namespace VID {
    enum VectorMesonType
    {
        JPSI,
        OMEGA,
        RHO,
        RHO_OMEGA_MIXING
    };
}

// Description:
//   Class to generate LambdaB -> Lambda(p pi) V(Vp Vm) decays
//   with V a vector meson such as J/psi (mu+mu-)
//                                 Rho (pi+pi-)
//                                 Omega (pi+pi-)
//                                 Rho-omega mixing (pi+pi-)
//
// DECAY : LambdaB -> Lambda + vector meson
//
// d(Sigma)
// -------- = 1 + A*B*cos(theta) + 2*A*Re(C*exp(i*phi))*sin(theta)
// d(Omega)
//
// with A (real)    : lambdaB  helicity asymmetry parameter
//      B (real)    : lambdaB polarisation
//      C (complex) : lambdaB density matrix element rho+-
//
// cf : O. Leitner, Z.J Ajaltouni, E. Conte,
//      PCCF RI 0601, ECT-05-15, LPNHE/2006-01, hep-ph/0602043

class EvtLambdaB2LambdaV : public EvtDecayProb {
  public:
    EvtLambdaB2LambdaV();

    EvtDecayBase* clone() override;

    std::string getName() override;
    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* lambdab ) override;

  private:
    //class name for report method
    std::string fname;

    //meson vector identity
    VID::VectorMesonType Vtype;

    //decay dynamics parameters
    double A;
    double B;
    EvtComplex C;

    //V mass generator method
    double getVMass( double MASS_LAMBDAB, double MASS_LAMBDA );

    //PDF generator method
    double BreitWignerRelPDF( double m, double _m0, double _g0 );
    double RhoOmegaMixingPDF( double m, double _mr, double _gr, double _mo,
                              double _go );
};

//*******************************************************************
//*                                                                 *
//*             Class EvtLambda2PPiForLambdaB2LambdaV               *
//*                                                                 *
//*******************************************************************
//
// DECAY : Lambda -> p + pi-
//
// d(Sigma)
// -------- = 1 + A*B*cos(theta) + 2*A*Re(D*exp(i*phi))*sin(theta)
// d(Omega)
//
// with A (real)    : lambda asymmetry parameter
//      B (real)    : lambda polarisation
//      C (real)    : lambdaB polarisation
//      D (complex) : lambda density matrix element rho+-
//
// cf : O. Leitner, Z.J Ajaltouni, E. Conte
//      PCCF RI 0601, ECT-05-15, LPNHE/2006-01, hep-ph/0602043

class EvtLambda2PPiForLambdaB2LambdaV : public EvtDecayProb {
  public:
    EvtLambda2PPiForLambdaB2LambdaV();
    EvtDecayBase* clone() override;

    std::string getName() override;
    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* lambda ) override;

  private:
    //class name for report method
    std::string fname;

    //meson vector identity
    VID::VectorMesonType Vtype;

    //decay dynamics parameters
    double A;
    double B;
    double C;
    EvtComplex D;
};

//*******************************************************************
//*                                                                 *
//*               Class EvtV2VpVmForLambdaB2LambdaV                 *
//*                                                                 *
//*******************************************************************
//
// DECAY : vector meson V -> Vp + Vm
//
// d(Sigma)
// -------- = (1-3A)*cos(theta)^2 + (1+A)   //leptonic decays
// d(Omega)
//
// d(Sigma)
// -------- = (3A-1)*cos(theta)^2 + (1-A)   //hadronic decays
// d(Omega)
//
// with A (real)    : V density matrix element indicating the
//                    probability to be longitudinally polarized
//
// cf : O. Leitner, Z.J Ajaltouni, E. Conte
//      PCCF RI 0601, ECT-05-15, LPNHE/2006-01, hep-ph/0602043

class EvtV2VpVmForLambdaB2LambdaV : public EvtDecayProb {
  public:
    EvtV2VpVmForLambdaB2LambdaV();

    EvtDecayBase* clone() override;

    std::string getName() override;
    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* V ) override;

  private:
    //class name for report method
    std::string fname;

    //meson vector identity
    VID::VectorMesonType Vtype;
    //decay dynamics parameters
    double A;
};

#endif
