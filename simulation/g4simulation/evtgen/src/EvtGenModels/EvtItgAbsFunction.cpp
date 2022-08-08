
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

#include "EvtGenModels/EvtItgAbsFunction.hh"

#include "EvtGenBase/EvtPatches.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}
#include "EvtGenBase/EvtReport.hh"

#include "assert.h"
using std::endl;

EvtItgAbsFunction::EvtItgAbsFunction( double lowerRange, double upperRange ) :
    _upperRange( upperRange ), _lowerRange( lowerRange )
{
}

double EvtItgAbsFunction::value( double x ) const
{
    if ( x >= _lowerRange && x <= _upperRange )
        return myFunction( x );
    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
        << "Error in EvtItgAbsFunction::value.  Given co-ordinate " << x
        << " is outside of allowed range [" << _lowerRange << ", "
        << _upperRange << "].  Returning 0.0" << endl;
    return 0.0;    // Never get here
}

double EvtItgAbsFunction::operator()( double x ) const
{
    return myFunction( x );
}
