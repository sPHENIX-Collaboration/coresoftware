
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

#ifndef EVTBTOSLLSCALARAMP_HH
#define EVTBTOSLLSCALARAMP_HH

#include "EvtGenModels/EvtbTosllAmp.hh"
class EvtParticle;
class EvtAmp;
class bTosllFF;

class EvtbTosllScalarAmp : public EvtbTosllAmp {
  public:
    //EvtbTosllScalarAmp(double c7, double c9, double c10):_c7(c7),_c9(c9),_c10(c10){}

    //Daughters are initialized and have been added to the parent.
    //No need to carry around the daughters seperately!
    void CalcAmp( EvtParticle* parent, EvtAmp& amp,
                  EvtbTosllFF* formFactors ) override;

  private:
    //double _c7,_c9,_c10;
};

#endif
