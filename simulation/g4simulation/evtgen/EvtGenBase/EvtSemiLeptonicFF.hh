
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

#ifndef EVTSEMILEPTONICFF_HH
#define EVTSEMILEPTONICFF_HH

#include "EvtGenBase/EvtId.hh"

class EvtSemiLeptonicFF {
  public:
    virtual ~EvtSemiLeptonicFF(){};

    virtual void getscalarff( EvtId parent, EvtId daught, double t, double mass,
                              double* fpf, double* fmf ) = 0;
    virtual void getvectorff( EvtId parent, EvtId daught, double t, double mass,
                              double* a1f, double* a2f, double* vf,
                              double* a0f ) = 0;
    virtual void gettensorff( EvtId parent, EvtId daught, double t, double mass,
                              double* a1f, double* a2f, double* vf,
                              double* a0f ) = 0;

    virtual void getbaryonff( EvtId parent, EvtId daught, double t,
                              double m_meson, double* f1v, double* f1a,
                              double* f2v, double* f2a ) = 0;

    virtual void getdiracff( EvtId parent, EvtId daught, double q2, double mass,
                             double* f1, double* f2, double* f3, double* g1,
                             double* g2, double* g3 ) = 0;

    virtual void getraritaff( EvtId parent, EvtId daught, double q2,
                              double mass, double* f1, double* f2, double* f3,
                              double* f4, double* g1, double* g2, double* g3,
                              double* g4 ) = 0;
};

#endif
