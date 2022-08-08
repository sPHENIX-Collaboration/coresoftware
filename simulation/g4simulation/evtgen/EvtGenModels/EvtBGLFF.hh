
/***********************************************************************
* Copyright 1998-2021 CERN for the benefit of the EvtGen authors       *
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

#ifndef EvtBGLFF_HH
#define EvtBGLFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

/** The class provides the form factors for semileptonic D and D* decays with full mass dependence
 */
class EvtBGLFF : public EvtSemiLeptonicFF {
  public:
    /** Default constructor */
    EvtBGLFF( double bglap_0, double bglap_1, double bglap_2, double bglap_3,
              double bgla0_0, double bgla0_1, double bgla0_2, double bgla0_3 );

    /** Default constructor */
    EvtBGLFF( double bgla_0, double bgla_1, double bglb_0, double bglb_1,
              double bglc_1, double bglc_2 );

    /** Returns scalar ffs */
    void getscalarff( EvtId parent, EvtId daught, double t, double mass,
                      double* fp, double* f0 ) override;

    /** Returns vector ffs */
    void getvectorff( EvtId parent, EvtId daught, double t, double mass,
                      double* a1f, double* a2f, double* vf,
                      double* a0f ) override;

    /** Returns tensor ffs */
    void gettensorff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    /** Returns baryon ffs */
    void getbaryonff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    /** Returns dirac ffs */
    void getdiracff( EvtId, EvtId, double, double, double*, double*, double*,
                     double*, double*, double* ) override;

    /** Returns tarita ffs */
    void getraritaff( EvtId, EvtId, double, double, double*, double*, double*,
                      double*, double*, double*, double*, double* ) override;

  private:
    /** B -> Dlnu:
      ai_n (i = p ---vector, 0 ---scalar; n = 0,1,2,3) are free coefficients of z expansion
      in dispersion relation parametrization from
      C.G.Boyd, B.Grinstein, R.F.Lebed, Phys. Rev. Lett. 74,4603(1995)

      Chosen the order of series N=3, i.e.
      a_0 + a_1 * z + a_2 * z^2 + a_3 * z^3

      Fitted values cited from
      R.Glattauer, etc. (Belle) Phys. Rev. D 93,032006 (2016).

      B -> D*lnu (l=e, mu):
      a_n, b_n (n = 0,1) and c_n (n = 0,1,2) are free coefficients of z expansion parametrization from
      C.G.Boyd, B.Grinstein and R.F.Lebed, Phys. Rev. D 56,6895(1997) &
      B.Grinstein, A.Kobach, Phys. Lett. B 771(2017)359-364

      For the expansion of form factors g and f, the order of series N=1, i.e.
      a_0 + a_1*z
      For the expansion of form factors F1, the order of series N=2, i.e.
      c_0 + c_1 * z + c_2 * z**2
      (g, f and F1 are the sub-terms of helicity amplitude)

      Fitted values are taken from a private discussion of Florian Bernlochner based on
      B.Grinstein and A.Kobach, Phys. Lett. B 771(2017)359-364

      It should not be used to generate D* with tau, due to the lack of fitted parameters in a0f amplitude.

  **/

    /** 0th-order z expansion coeffieient for vector form factor: f_+  */
    double m_ap_0{ 0 };

    /** 1st-order z expansion coeffieient for vector form factor: f_+  */
    double m_ap_1{ 0 };

    /** 2nd-order z expansion coeffieient for vector form factor: f_+  */
    double m_ap_2{ 0 };

    /** 3rd-order z expansion coeffieient for vector form factor: f_+  */
    double m_ap_3{ 0 };

    /** 0th-order z expansion coeffieient for scalar form factor f_0   */
    double m_a0_0{ 0 };

    /** 1st-order z expansion coeffieient for scalar form factor f_0   */
    double m_a0_1{ 0 };

    /** 2nd-order z expansion coeffieient for scalar form factor f_0   */
    double m_a0_2{ 0 };

    /** 3rd-order z expansion coeffieient for scalar form factor f_0   */
    double m_a0_3{ 0 };

    /** B->D*lnu z expansion coefficients  */

    /** 0th-order z expansion coefficient for form factor g   */
    double m_a_0{ 0 };

    /** 1st-order z expansion coefficient for form factor g   */
    double m_a_1{ 0 };

    /** 0th-order z expansion coefficient for form factor f   */
    double m_b_0{ 0 };

    /** 1st-order z expansion coefficient for form factor f   */
    double m_b_1{ 0 };

    /** 1st-order z expansion coefficient for form factor F1   */
    double m_c_1{ 0 };

    /** 2nd-order z expansion coefficient for form factor F1   */
    double m_c_2{ 0 };
};
#endif
