
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

#ifndef EVTKKLAMBDACFF_HH
#define EVTKKLAMBDACFF_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtKKLambdaCFF : public EvtSemiLeptonicFF {
  public:
    EvtKKLambdaCFF( int numarg, double* arglist );

    void getscalarff( EvtId parent, EvtId daught, double t, double mass,
                      double* f0p, double* f0m ) override;

    void getvectorff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    void gettensorff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    void getbaryonff( EvtId parent, EvtId daught, double t, double mass,
                      double* f1v, double* f1a, double* f2v,
                      double* f2a ) override;

    void getdiracff( EvtId, EvtId, double, double, double*, double*, double*,
                     double*, double*, double* ) override;

    void getraritaff( EvtId, EvtId, double, double, double*, double*, double*,
                      double*, double*, double*, double*, double* ) override;

  private:
    int _nargs;
    double _args[2];
};

#endif
