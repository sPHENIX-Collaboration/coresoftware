
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

#ifndef EVTDFUNCTIONSINGLE_HH
#define EVTDFUNCTIONSINGLE_HH

// Description:Evaluation of one Wigner d-Functions

class EvtdFunctionSingle {
  public:
    EvtdFunctionSingle();
    ~EvtdFunctionSingle();

    void init( int j, int m1, int m2 );

    double d( int j, int m1, int m2, double theta );

  private:
    int fact( int n );

    int _j;
    int _m1;
    int _m2;

    double* _coef;

    int _kmin;
    int _kmax;
};

#endif
