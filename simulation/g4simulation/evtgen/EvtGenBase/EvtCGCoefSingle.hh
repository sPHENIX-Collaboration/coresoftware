
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

#ifndef EVTCGCOEFSINGLE_HH
#define EVTCGCOEFSINGLE_HH

#include <vector>

class EvtCGCoefSingle final {
  public:
    EvtCGCoefSingle( int j1, int j2 ) { init( j1, j2 ); }

    double coef( int J, int M, int j1, int j2, int m1, int m2 );

  private:
    void init( int j1, int j2 );
    double& cg( int J, int M, int m1, int m2 );

    int _j1;
    int _j2;

    int _Jmax;
    int _Jmin;

    std::vector<std::vector<std::vector<double>>> _table;
};

#endif
