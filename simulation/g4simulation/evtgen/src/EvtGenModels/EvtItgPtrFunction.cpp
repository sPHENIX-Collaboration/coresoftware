
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

#include "EvtGenModels/EvtItgPtrFunction.hh"

#include "EvtGenBase/EvtPatches.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//----------------
// Constructors --
//----------------
EvtItgPtrFunction::EvtItgPtrFunction(
    double ( *theFunction )( double, const std::vector<double>& ),
    double lowerRange, double upperRange, const std::vector<double>& coeffs1 ) :
    EvtItgAbsFunction( lowerRange, upperRange ),
    _myFunction( theFunction ),
    _coeffs1( coeffs1 )
{
}

double EvtItgPtrFunction::myFunction( double x ) const
{
    return _myFunction( x, _coeffs1 );
}

void EvtItgPtrFunction::setCoeff( int vect, int which, double value )
{
    if ( vect == 1 )
        _coeffs1[which] = value;
}

double EvtItgPtrFunction::getCoeff( int vect, int which )
{
    if ( vect == 1 )
        return _coeffs1[which];
    else {
        return 0;
    }
}
