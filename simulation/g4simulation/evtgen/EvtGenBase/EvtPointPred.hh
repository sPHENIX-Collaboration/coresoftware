
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

#ifndef EVT_POINT_PRED_HH
#define EVT_POINT_PRED_HH

// Predicate testing validity of a point. The point class must provide
// bool isValid() method

template <class Point>
class EvtPointPred {
  public:
    typedef Point argument_type;
    typedef bool result_type;

    EvtPointPred() {}
    EvtPointPred( const EvtPointPred& ) {}
    ~EvtPointPred() {}

    result_type operator()( argument_type x ) { return x.isValid(); }
};
