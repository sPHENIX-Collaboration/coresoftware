
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

#ifndef EVTBTOKPIPICP_HH
#define EVTBTOKPIPICP_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtBTo3hCP.hh"

class EvtParticle;

// Description: Routine to decay B->K pi pi
//              and has CP violation.
//       --- This is the routine to be called by the Main generator
//          to get the decay of B0    -->-- K+ pi- pi0
//          The decay proceeeds through three channels:
//          a) B0 -->-- K*+ pi-  ; K*+    -->-- K+ pi0
//          b)          K*0 pi0  ; K*0bar -->-- K+ pi-
//          c)          K-  rho+ ; rho+   -->-- pi+ pi0
//         It provides at the same time the CP conjugate decay
//                              B0bar -->-- K- pi+ pi0

class EvtBToKpipiCP : public EvtDecayAmp {
  public:
    EvtBToKpipiCP() {}

    std::string getName() override;
    EvtBToKpipiCP* clone() override;

    void init() override;
    void decay( EvtParticle* p ) override;

  private:
    EvtBTo3hCP generator;
};

#endif
