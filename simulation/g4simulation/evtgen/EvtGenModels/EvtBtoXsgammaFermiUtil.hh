
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

#ifndef EVTBTOXSGAMMAFERMIUTIL_HH
#define EVTBTOXSGAMMAFERMIUTIL_HH

#include <vector>

// Description:
//       Class to hold various fermi functions and their helper functions. The
//       fermi functions are used in EvtBtoXsgammaKagan.

class EvtBtoXsgammaFermiUtil final {
    //--------------------
    // Instance Members --
    //--------------------

  public:
    //Exponential function
    static double FermiExpFunc( double var, const std::vector<double>& coeffs );

    //Gaussian function and its helper functions
    static double FermiGaussFunc( double, std::vector<double> const& coeffs );
    static double FermiGaussFuncRoot( double, double, double,
                                      std::vector<double>& coeffs );
    static double FermiGaussRootFcnA( double, const std::vector<double>& coeffs1,
                                      const std::vector<double>& coeffs2 );
    static double FermiGaussRootFcnB( double, const std::vector<double>& coeffs1,
                                      const std::vector<double>& coeffs2 );
    static double Gamma( double, const std::vector<double>& coeffs );

    //Roman function and its helper functions
    static double BesselI1( double );
    static double BesselK1( double );
    static double FermiRomanFuncRoot( double, double );
    static double FermiRomanRootFcnA( double );
    static double FermiRomanFunc( double, std::vector<double> const& coeffs );
};

#endif
