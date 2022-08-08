
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

#ifndef EVTITGFUNCTION_HH
#define EVTITGFUNCTION_HH

#include "EvtGenModels/EvtItgAbsFunction.hh"

/**
 *  Generic function where the pointer to the function is available.
 *
 *  The function is taken as type pointer to function returning double and
 *  taking a double (the abscissa) and a const RWTValVector<double> reference
 *  (the parameter values of the function) as arguments.
 */

class EvtItgFunction : public EvtItgAbsFunction {
  public:
    // Constructors
    EvtItgFunction( double ( *theFunction )( double ), double lowerRange,
                    double upperRange );

    void setCoeff( int, int, double ) override{};
    double getCoeff( int, int ) override { return 0.0; };

  protected:
    // Helper functions

    double myFunction( double x ) const override;

  private:
    // Data members
    double ( *_myFunction )( double x );
};

#endif    // EvtITGFUNCTION_HH
