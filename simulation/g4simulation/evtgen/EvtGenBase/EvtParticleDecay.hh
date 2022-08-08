
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

#ifndef EVTPARTICLEDECAY_HH
#define EVTPARTICLEDECAY_HH

#include "EvtGenBase/EvtDecayBase.hh"

class EvtParticleDecay {
  public:
    EvtParticleDecay()
    {
        _decay = 0;
        _brfrsum = 0.0;
        _massmin = 0.0;
    }

    ~EvtParticleDecay()
    {
        if ( _decay != 0 )
            delete _decay;
    }

    void chargeConj( EvtParticleDecay* decay );

    void setDecayModel( EvtDecayBase* decay ) { _decay = decay; }
    EvtDecayBase* getDecayModel() { return _decay; }
    double getBrfrSum() { return _brfrsum; }
    void setBrfrSum( double brfrsum ) { _brfrsum = brfrsum; }
    double getMassMin() { return _massmin; }
    void setMassMin( double massmin ) { _massmin = massmin; }

    void printSummary();

  private:
    EvtDecayBase* _decay;

    double _brfrsum;
    double _massmin;
};

#endif
