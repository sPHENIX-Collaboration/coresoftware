
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

#ifndef EVTPARTICLEDECAYLIST_HH
#define EVTPARTICLEDECAYLIST_HH

#include "EvtGenBase/EvtParticleDecay.hh"

typedef EvtParticleDecay* EvtParticleDecayPtr;

class EvtParticleDecayList {
  public:
    EvtParticleDecayList()
    {
        _decaylist = 0;
        _nmode = 0;
        _rawbrfrsum = 0;
    }

    EvtParticleDecayList( const EvtParticleDecayList& o );

    ~EvtParticleDecayList();

    EvtParticleDecayList& operator=( const EvtParticleDecayList& o );

    int getNMode() const { return _nmode; }

    void setNMode( int nmode );

    EvtDecayBase* getDecayModel( EvtParticle* p );
    EvtDecayBase* getDecayModel( int imode );

    EvtParticleDecay& getDecay( int nchannel ) const;

    double getRawBrfrSum() { return _rawbrfrsum; }
    void setRawBrfrSum( double rawbrfrsum ) { _rawbrfrsum = rawbrfrsum; }

    void makeChargeConj( EvtParticleDecayList* conjDecayList );

    void removeDecay();

    void alocateDecay( int nmode )
    {
        _decaylist = new EvtParticleDecayPtr[nmode];
    }

    void removeMode( EvtDecayBase* decay );

    void addMode( EvtDecayBase* decay, double brfr, double massmin );
    void finalize();

    void printSummary();

    bool isJetSet() const;

  private:
    EvtParticleDecayPtr* _decaylist;

    double _rawbrfrsum;
    int _nmode;
};

#endif
