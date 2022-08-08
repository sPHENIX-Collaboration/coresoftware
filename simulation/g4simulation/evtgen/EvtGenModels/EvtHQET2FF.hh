
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

#ifndef EVTHQET2FF_HH
#define EVTHQET2FF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class EvtHQET2FF : public EvtSemiLeptonicFF {
  public:
    EvtHQET2FF( double hqetrho2, double hqetha1_1, double hqetr1_1,
                double hqetr2_1 );
    EvtHQET2FF( double hqetrho2, double hqetv1_1 );
    EvtHQET2FF( double hqetrho2, double hqetha1_1, double hqetr1_1,
                double hqetr2_1, double hqetr0_1 );
    EvtHQET2FF( double hqetrho2, double hqetv1_1, double indelta );
    void getvectorff( EvtId parent, EvtId daught, double t, double mass,
                      double* a1f, double* a2f, double* vf,
                      double* a0f ) override;

    void getscalarff( EvtId parent, EvtId daught, double t, double mass,
                      double* f0p, double* f0m ) override;

    void gettensorff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    void getbaryonff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    void getdiracff( EvtId, EvtId, double, double, double*, double*, double*,
                     double*, double*, double* ) override;

    void getraritaff( EvtId, EvtId, double, double, double*, double*, double*,
                      double*, double*, double*, double*, double* ) override;

  private:
    double r1_1;
    double rho2;
    double r2_1;
    double ha1_1;
    double v1_1;
    double r0_1;
    double delta;
    bool extended;
};

#endif
