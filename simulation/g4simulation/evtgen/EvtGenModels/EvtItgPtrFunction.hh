
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

#ifndef EVTITGPTRFUNCTION_HH
#define EVTITGPTRFUNCTION_HH

#include "EvtGenModels/EvtItgAbsFunction.hh"

#include <vector>

// Description:
//      Class describing a function with one vector of coefficients. (Stolen and
//      modified from the BaBar IntegrationUtils package - author: Phil Strother).

class EvtItgPtrFunction : public EvtItgAbsFunction {
  public:
    EvtItgPtrFunction( double ( *theFunction )( double,
                                                const std::vector<double>& ),
                       double lowerRange, double upperRange,
                       const std::vector<double>& coeffs1 );

    void setCoeff( int, int, double ) override;
    double getCoeff( int, int ) override;

  protected:
    double myFunction( double x ) const override;

  private:
    // Data members
    double ( *_myFunction )( double x, const std::vector<double>& coeffs1 );

    std::vector<double> _coeffs1;
};

#endif    // EVTITGPTRFUNCTION_HH
