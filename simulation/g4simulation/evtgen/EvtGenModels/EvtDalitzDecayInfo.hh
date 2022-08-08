
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

#ifndef EVTDALITZDECAYINFO_HH
#define EVTDALITZDECAYINFO_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtId.hh"

#include <vector>

// Description: Model to describe a generic dalitz decay

class EvtDalitzDecayInfo final {
  public:
    EvtDalitzDecayInfo( EvtId d1, EvtId d2, EvtId d3 ) :
        _d1( d1 ), _d2( d2 ), _d3( d3 ), _probMax( 0. )
    {
    }

    void addResonance( EvtComplex amp, EvtDalitzReso res )
    {
        _resonances.push_back( std::pair<EvtComplex, EvtDalitzReso>( amp, res ) );
    }
    void addResonance( std::pair<EvtComplex, EvtDalitzReso> res )
    {
        _resonances.push_back( res );
    }
    void setProbMax( double probMax ) { _probMax = probMax; }

    const std::vector<std::pair<EvtComplex, EvtDalitzReso>>& getResonances() const
    {
        return _resonances;
    }
    double getProbMax() const { return _probMax; }

    inline const EvtId& daughter1() const { return _d1; }
    inline const EvtId& daughter2() const { return _d2; }
    inline const EvtId& daughter3() const { return _d3; }

  private:
    EvtId _d1, _d2, _d3;
    std::vector<std::pair<EvtComplex, EvtDalitzReso>> _resonances;
    double _probMax;
};

#endif
