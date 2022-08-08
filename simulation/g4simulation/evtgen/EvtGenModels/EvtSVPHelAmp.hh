
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

#ifndef EVTSVPHELAMP_HH
#define EVTSVPHELAMP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

// Description: Routine to decay scalar -> vector + photon
//              by specifying the helicity amplitudes
//
//Class to handle decays of the form SCALAR -> VECTOR PHOTON
//where the helicity amplitudes must be specified. The
//first and third arguments are the magnitudes of the H+
//and H- helicity amplitudes respectively. The second and
//fourth arguements are the phases.
//Calls EvtSVPHel.

class EvtParticle;
class EvtAmp;
class EvtId;

class EvtSVPHelAmp : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;

    static void SVPHel( EvtParticle* parent, EvtAmp& amp, EvtId n_v1,
                        EvtId n_ph, const EvtComplex& hp, const EvtComplex& hm );
};

#endif
