
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

#ifndef EVTBS2LLGAMMAFFMNT_HH
#define EVTBS2LLGAMMAFFMNT_HH

#include "EvtGenModels/Evtbs2llGammaFF.hh"

class EvtId;

class EvtComplex;

// Description: Form factors for rare radiative leptonic B-decays
//              according to the papers:
//              1) F.Kruger, D.Melikhov,  Phys. Rev. D67, 034002, 2003.
//              2) D.Melikhov, N.Nikitin, Phys. Rev. D70, 114028, 2004.
//              3) I.Balakireva, D.Melikhov, N.Nikitin, D.Tlisov,
//                                           e-Print: arXiv:0911.0605 [hep-ph].
// Description: Form-factors for B_q -> gamma transitions,  q={d,s}
//              according to the papers:
//              1) F.Kruger, D.Melikhov,  Phys. Rev. D67, 034002, 2003.
//              2) D.Melikhov, N.Nikitin, Phys. Rev. D70, 114028, 2004.
//              3) I.Balakireva, D.Melikhov, N.Nikitin, D.Tlisov,
//                                        Phys. Rev. D81, 054024, 2010.

class Evtbs2llGammaFFMNT : public Evtbs2llGammaFF {
  public:
    Evtbs2llGammaFFMNT();

    void getPhotonFF( int decay_id, double fb, EvtId parent, double q2,
                      double M1, double mb, double mq, EvtComplex c7gam,
                      EvtComplex a1, EvtComplex lambda_qu, EvtComplex lambda_qc,
                      EvtComplex& Fv, EvtComplex& Fa, EvtComplex& Ftv,
                      EvtComplex& Fta ) override;

    double getQuarkMass( int i ) override;
};

#endif
