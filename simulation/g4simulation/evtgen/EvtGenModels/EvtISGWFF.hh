
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

#ifndef EVTISGWFF_HH
#define EVTISGWFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

// Description:Form factor routines specific to EvtISGW

class EvtISGWFF : public EvtSemiLeptonicFF {
    void getscalarff( EvtId parent, EvtId daught, double t, double mass,
                      double* fpf, double* f0f ) override;
    void getvectorff( EvtId parent, EvtId daught, double t, double mass,
                      double* a1f, double* a2f, double* vf,
                      double* a0f ) override;
    void gettensorff( EvtId parent, EvtId daught, double t, double mass,
                      double* hf, double* kf, double* bpf, double* bmf ) override;

    void getbaryonff( EvtId, EvtId, double, double, double*, double*, double*,
                      double* ) override;

    void getdiracff( EvtId, EvtId, double, double, double*, double*, double*,
                     double*, double*, double* ) override;

    void getraritaff( EvtId, EvtId, double, double, double*, double*, double*,
                      double*, double*, double*, double*, double* ) override;

    // getscalarff, getvectorff, and gettensorff call the
    // correct isgw form factor routine which computes
    // form factors according to the ISGW paper.

    void EvtISGW1FF3S1( EvtId parent, EvtId daught, double t, double mass,
                        double* ff, double* gf, double* apf, double* amf );
    void EvtISGW1FF23S1( EvtId parent, EvtId daught, double t, double mass,
                         double* fpf, double* gpf, double* app, double* apm );
    void EvtISGW1FF3P1( EvtId parent, EvtId daught, double t, double mass,
                        double* lf, double* qf, double* cpf, double* cmf );
    void EvtISGW1FF3P0( EvtId parent, EvtId daught, double t, double mass,
                        double* upf, double* umf );
    void EvtISGW1FF1S0( EvtId parent, EvtId daught, double t, double mass,
                        double* fpf, double* fmf );
    void EvtISGW1FF21S0( EvtId parent, EvtId daught, double t, double mass,
                         double* fppf, double* fpmf );
    void EvtISGW1FF3P2( EvtId parent, EvtId daught, double t, double mass,
                        double* h, double* k, double* bp, double* bm );
    void EvtISGW1FF1P1( EvtId parent, EvtId daught, double t, double mass,
                        double* vf, double* rf, double* spf, double* smf );
};

#endif
