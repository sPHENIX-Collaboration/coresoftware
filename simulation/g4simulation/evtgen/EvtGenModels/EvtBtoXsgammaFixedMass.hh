
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

#ifndef EVTBTOXSGAMMAFIXEDMASS_HH
#define EVTBTOXSGAMMAFIXEDMASS_HH

#include "EvtGenModels/EvtBtoXsgammaAbsModel.hh"

// Description:
//       Implimentation of a fixed hadronic mass to measure spectrum

class EvtBtoXsgammaFixedMass : public EvtBtoXsgammaAbsModel {
  public:
    void init( int, double* ) override;

    double GetMass( int code ) override;

  private:
    //Input parameters
    double _mH;
};

#endif
