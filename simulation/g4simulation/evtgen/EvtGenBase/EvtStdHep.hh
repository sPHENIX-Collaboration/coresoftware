
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

#ifndef EVTSTDHEP_HH
#define EVTSTDHEP_HH

#include "EvtGenBase/EvtVector4R.hh"

#include <iosfwd>

const int EVTSTDHEPLENGTH = 1000;

class EvtStdHep {
  public:
    EvtStdHep() {}
    ~EvtStdHep() {}

    void init();

    int getFirstMother( int i ) { return _prntfirst[i]; }
    int getLastMother( int i ) { return _prntlast[i]; }
    int getFirstDaughter( int i ) { return _daugfirst[i]; }
    int getLastDaughter( int i ) { return _dauglast[i]; }

    int getStdHepID( int i ) { return _id[i]; }
    int getIStat( int i ) { return _istat[i]; }

    EvtVector4R getP4( int i ) { return _p4[i]; }
    EvtVector4R getX4( int i ) { return _x[i]; }

    void translate( EvtVector4R d );

    int getNPart();
    void createParticle( EvtVector4R p4, EvtVector4R x, int prntfirst,
                         int prntlast, int id );

    friend std::ostream& operator<<( std::ostream& s, const EvtStdHep& stdhep );

  private:
    int _npart;
    EvtVector4R _p4[EVTSTDHEPLENGTH];
    EvtVector4R _x[EVTSTDHEPLENGTH];
    int _prntfirst[EVTSTDHEPLENGTH];
    int _prntlast[EVTSTDHEPLENGTH];
    int _daugfirst[EVTSTDHEPLENGTH];
    int _dauglast[EVTSTDHEPLENGTH];
    int _id[EVTSTDHEPLENGTH];
    int _istat[EVTSTDHEPLENGTH];
};

#endif
