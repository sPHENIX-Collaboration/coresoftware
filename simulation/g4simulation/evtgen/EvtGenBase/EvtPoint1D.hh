
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

#ifndef EVT_POINT_1D_HH
#define EVT_POINT_1D_HH

// Point on a finite 1-D interval. isValid shows whether for a given specification,
// the coordinate _value is inside the interval defined by _min, _max.

class EvtPoint1D final {
  public:
    EvtPoint1D();
    EvtPoint1D( double value );
    EvtPoint1D( double min, double max, double value );

    bool isValid() const { return _valid; }

    double value() const { return _value; }

    void print() const;

  private:
    double _min;    // interval minimum
    double _max;    // interval maximum
    double _value;
    bool _valid;    // valid point inside the interval?
};

#endif
