
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

#ifndef EVTSVVHELAMP_HH
#define EVTSVVHELAMP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

//Class to handle decays of the form SCALAR -> VECTOR VECTOR
//according the the helicity amplitudes specified by the
//user.  There are 6 arguements, orders as amplitude then
//phase for H+, H0, and H-, in that order.

class EvtAmp;
class EvtParticle;
class EvtId;

class EvtSVVHelAmp : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;

    static void SVVHel( EvtParticle* parent, EvtAmp& amp, EvtId n_v1,
                        EvtId n_v2, const EvtComplex& hp, const EvtComplex& h0,
                        const EvtComplex& hm );

    std::string getParamName( int i ) override;
    std::string getParamDefault( int i ) override;
};

#endif
