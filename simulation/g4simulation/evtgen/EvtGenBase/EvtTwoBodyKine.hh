
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

#ifndef EVT_TWO_BODY_KINE_HH
#define EVT_TWO_BODY_KINE_HH

#include <iostream>

// Descriptions of the kinematics of a two-body decay.

class EvtTwoBodyKine {
  public:
    enum Index
    {
        A,
        B,
        AB
    };

    EvtTwoBodyKine();
    EvtTwoBodyKine( double mA, double mB, double mAB );

    // Accessors

    inline double mA() const { return _mA; }
    inline double mB() const { return _mB; }
    inline double mAB() const { return _mAB; }
    double m( Index i ) const;

    // Momentum of the other two particles in the
    // rest-frame of particle i.

    double p( Index i = AB ) const;

    // Energy of particle i in the rest frame of particle j

    double e( Index i, Index j ) const;

    void print( std::ostream& os ) const;

  private:
    double _mA;
    double _mB;
    double _mAB;
};

std::ostream& operator<<( std::ostream& os, const EvtTwoBodyKine& p );

#endif
