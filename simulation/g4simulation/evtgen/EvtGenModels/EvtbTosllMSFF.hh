
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

#ifndef EVTBTOSLLMSFF_HH
#define EVTBTOSLLMSFF_HH

#include "EvtGenModels/EvtbTosllFFNew.hh"

class EvtId;

// Description: Form factors for B^0_d -> K^* l^+ l^- transition according
//              to the paper: D.Melikhov, B.Stech, PRD62, 014006 (2000).
//  N.Nikitin (nnikit@mail.cern.ch)  March 13, 2008   Module created
//  N.Nikitin (nnikit@mail.cern.ch)  March 27, 2008   add \bar B -> \bar (K,K^*) transition ff
//  N.Nikitin (nnikit@mail.cern.ch)  April 26, 2008   add \bar Bs -> phi transition ff
//  N.Nikitin (nnikit@mail.cern.ch)  April 26, 2008   add \bar Bs -> K*  transition ff
//  N.Nikitin (nnikit@mail.cern.ch)  April 27, 2008   add \bar B -> \bar rho transition ff
//  N.Nikitin (nnikit@mail.cern.ch)  Nvmbr 04, 2011   add \bar B -> omega transition ff
//  N.Nikitin (nnikit@mail.cern.ch)  Dec   16, 2011   add \bar B -> \bar K_1(1270) transition ff (from H.Hatanaka and Kwei-Chou Yang, PRD78, 074007 (2008))
//  N.Nikitin (nnikit@mail.cern.ch)  Dec   16, 2011   add \bar B -> \bar K_1(1400) transition ff (from H.Hatanaka and Kwei-Chou Yang, PRD78, 074007 (2008))
//  N.Nikitin (nnikit@mail.cern.ch)  May   11, 2012   add \bar Bs -> f_0(980) transition ff with NLO corrections from Table II (see P.Colangelo et al., PRD81, 074001 (2010))

class EvtbTosllMSFF : public EvtbTosllFFNew {
  public:
    EvtbTosllMSFF();

    double equation9_10( double ff0, double M2, double q2, double sigma1,
                         double sigma2, int eq_num );

    void getScalarFF( EvtId parent, EvtId daught, double t, double& fp,
                      double& f0, double& ft ) override;

    void getVectorFF( EvtId parent, EvtId daught, double t, double& a1,
                      double& a2, double& a0, double& v, double& t1, double& t2,
                      double& t3 ) override;

    double getQuarkMass( int i ) override;

  private:
};

#endif
