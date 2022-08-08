
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

#ifndef EVTDTOKPIENU_HH
#define EVTDTOKPIENU_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <array>

class EvtParticle;

class EvtDToKpienu : public EvtDecayProb {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

  private:
    void KinVGen( const EvtVector4R& vp4_K, const EvtVector4R& vp4_Pi,
                  const EvtVector4R& vp4_Lep, const EvtVector4R& vp4_Nu,
                  const int charm, double& m2, double& q2, double& cosV,
                  double& cosL, double& chi ) const;
    double calPDF( const double m2, const double q2, const double cosV,
                   const double cosL, const double chi ) const;
    void ResonanceP( const double m, const double q, const double mV,
                     const double mA, const double V_0, const double A1_0,
                     const double A2_0, const double m0, const double width0,
                     const double rBW, double& amplitude, double& delta,
                     EvtComplex& F11, EvtComplex& F21, EvtComplex& F31 ) const;
    void NRS( const double m, const double q, const double rS, const double rS1,
              const double a_delta, const double b_delta, const double mA,
              const double m0, const double width0, double& amplitude,
              double& delta, EvtComplex& F10 ) const;
    void ResonanceD( const double m, const double q, const double mV,
                     const double mA, const double TV_0, const double T1_0,
                     const double T2_0, const double m0, const double width0,
                     const double rBW, double& amplitude, double& delta,
                     EvtComplex& F12, EvtComplex& F22, EvtComplex& F32 ) const;
    double getPStar( const double m, const double m1, const double m2 ) const;
    double getF1( const double m, const double m0, const double m_c1,
                  const double m_c2, const double rBW ) const;
    double getF2( const double m, const double m0, const double m_c1,
                  const double m_c2, const double rBW ) const;
    double getWidth0( const double m, const double m0, const double m_c1,
                      const double m_c2, const double width0 ) const;
    double getWidth1( const double m, const double m0, const double m_c1,
                      const double m_c2, const double width0,
                      const double rBW ) const;
    double getWidth2( const double m, const double m0, const double m_c1,
                      const double m_c2, const double width0,
                      const double rBW ) const;
    EvtComplex getCoef( const double rho, const double phi ) const;

    int nAmps;
    std::array<int, 5> type;

    double rS;
    double rS1;
    double a_delta;
    double b_delta;
    double m0_1430_S;
    double width0_1430_S;

    double mV;
    double mA;
    double V_0;
    double A1_0;
    double A2_0;
    double m0;
    double width0;
    double rBW;
    double rho;
    double phi;
    double m0_1410;
    double width0_1410;
    double rho_1410;
    double phi_1410;
    double TV_0;
    double T1_0;
    double T2_0;
    double m0_1430;
    double width0_1430;
    double rho_1430;
    double phi_1430;

    double mD;
    double mPi;
    double mK;
    double Pi;
    double root2;
    double root2d3;
    double root1d2;
    double root3d2;
};

#endif
