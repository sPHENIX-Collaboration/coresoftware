
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

#ifndef EVTBS2LLGAMMAFF_HH
#define EVTBS2LLGAMMAFF_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtId.hh"

// Description: This is the NEW base class for form factors in b->sll transitions.

class Evtbs2llGammaFF {
  public:
    virtual ~Evtbs2llGammaFF(){};

    virtual void getPhotonFF( int /*decay_id*/, double /*fb*/, EvtId /*parent*/,
                              double /*q2*/, double /*M1*/, double /*mb*/,
                              double /*mq*/, EvtComplex /*c7gam*/,
                              EvtComplex /*a1*/, EvtComplex /*lambda_qu */,
                              EvtComplex /*lambda_qc*/, EvtComplex& /*Fv*/,
                              EvtComplex& /*Fa*/, EvtComplex& /*Ftv*/,
                              EvtComplex& /*Fta*/ )
    {
        return;
    };

    virtual double getQuarkMass( int /*i*/ ) { return 0.0; };
};

#endif
