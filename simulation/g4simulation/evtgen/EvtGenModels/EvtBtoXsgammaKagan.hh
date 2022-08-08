
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

#ifndef EVTBTOXSGAMMAKAGAN_HH
#define EVTBTOXSGAMMAKAGAN_HH

#include "EvtGenModels/EvtBtoXsgammaAbsModel.hh"

#include <vector>

// Description:
//       Implimentation of the Kagan-Neubert model for non-resonant
//       B->Xs,gamma decays.
// Description:
//       Routine to perform two-body non-resonant B->Xs,gamma decays.
//       The X_s mass spectrum generated is based on the Kagan-Neubert model.
//       See hep-ph/9805303 for the model details and input parameters.
//
//       The input parameters are 1:fermi_model, 2:mB, 3:mb, 4:mu, 5:lam1,
//       6:delta, 7:z, 8:nIntervalS, 9:nIntervalmH. Choosing fermi_model=1
//       uses an exponential shape function, fermi_model=2 uses a gaussian
//       shape function and fermi_model=3 a roman shape function. The complete mass
//       spectrum for a given set of input parameters is calculated from
//       scratch in bins of nIntervalmH. The s22, s27 and s28 coefficients are calculated
//       in bins of nIntervalS. As the program includes lots of integration, the
//       theoretical hadronic mass spectra is computed for the first time
//       the init method is called. Then, all the other times (eg if we want to decay a B0
//       as well as an anti-B0) the vector mass info stored the first time is used again.

class EvtBtoXsgammaKagan : public EvtBtoXsgammaAbsModel {
  public:
    void init( int, double* ) override;

    void computeHadronicMass( int, double* );

    void getDefaultHadronicMass();

    double GetMass( int code ) override;

    double CalcAlphaS( double );

    void CalcWilsonCoeffs();
    void CalcDelta();
    double Fz( double );

  private:
    //Input parameters
    double _mb;
    double _mB;
    double _delta;
    double _nIntervalS;
    double _nIntervalmH;
    double _lambdabar;
    double _lam1;
    double _mHmin;
    double _mHmax;
    //Other parameters
    double _r7;
    double _gam77;
    double _gam27;
    double _gam87;
    double _beta0;
    double _beta1;
    double _alphasmZ;
    double _mZ;
    double _z;
    double _fz;
    double _lam2;
    double _kappabar;
    double _rer2;
    double _rer8;
    double _kSLemmu;
    double _mW;
    double _mt;
    double _ms;
    double _mu;

    double _c2mu;
    double _c70mu;
    double _c80mu;
    double _c71mu;
    double _c7emmu;

    double _cDeltatot;

    double _alpha;
    double _alphasmW;
    double _alphasmt;
    double _alphasmu;
    double _alphasmubar;
    double _etamu;

    std::vector<double> _mHVect;

    static double ReG( double );
    static double ImG( double );
    static double s77( double );
    static double s88( double, double, double );
    static double s78( double );
    static double s22Func( double var, const std::vector<double>& coeffs );
    static double s27Func( double var, const std::vector<double>& coeffs );

    static double Delta( double, double );
    static double DeltaFermiFunc( double, const std::vector<double>& coeffs1,
                                  const std::vector<double>& coeffs2,
                                  const std::vector<double>& coeffs3 );
    static double s77FermiFunc( double, const std::vector<double>& coeffs1,
                                const std::vector<double>& coeffs2 );
    static double s88FermiFunc( double, const std::vector<double>& coeffs1,
                                const std::vector<double>& coeffs2,
                                const std::vector<double>& coeffs3 );
    static double s78FermiFunc( double, const std::vector<double>& coeffs1,
                                const std::vector<double>& coeffs2 );
    static double s22FermiFunc( double, std::vector<double>& coeffs );
    static double s27FermiFunc( double, std::vector<double>& coeffs );
    static double s28FermiFunc( double, std::vector<double>& coeffs );
    static double GetArrayVal( double, double, double, double,
                               std::vector<double> );
    static double sFermiFunc( double, const std::vector<double>& coeffs1,
                              const std::vector<double>& coeffs2,
                              const std::vector<double>& coeffs3,
                              const std::vector<double>& coeffs4 );
    static double FermiFunc( double, const std::vector<double>& coeffs );
    static double diLogFunc( double );
    static double diLogMathematica( double );
    std::vector<double> massHad, brHad;
    static double intervalMH;
    static bool bbprod;
};

#endif
