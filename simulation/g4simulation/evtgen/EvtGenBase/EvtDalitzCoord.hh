
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

#ifndef EVT_DALITZ_COORD_HH
#define EVT_DALITZ_COORD_HH

#include "EvtGenBase/EvtCyclic3.hh"

#include <iostream>

// Two dimensional coordinate of a point in a Dalitz plot

class EvtDalitzCoord final {
  public:
    // ctor, dtor

    EvtDalitzCoord();
    EvtDalitzCoord( EvtCyclic3::Pair i1, double q1, EvtCyclic3::Pair i2,
                    double q2 );
    EvtDalitzCoord( const EvtDalitzCoord& other );

    inline EvtCyclic3::Pair pair1() const { return _i1; }
    inline EvtCyclic3::Pair pair2() const { return _i2; }
    inline double q1() const { return _q1; }
    inline double q2() const { return _q2; }

    // It's nice to have an equality operator for
    // a coordinate. However, beware effects of numerical precision

    bool operator==( const EvtDalitzCoord& ) const;

    void print( std::ostream& ) const;

  private:
    // Two coordinates define the point

    EvtCyclic3::Pair _i1;
    EvtCyclic3::Pair _i2;

    double _q1;
    double _q2;
};

std::ostream& operator<<( std::ostream&, const EvtDalitzCoord& );

#endif
