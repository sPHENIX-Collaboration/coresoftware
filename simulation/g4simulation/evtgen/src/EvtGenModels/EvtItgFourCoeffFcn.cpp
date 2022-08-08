
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

#include "EvtGenModels/EvtItgFourCoeffFcn.hh"

#include "EvtGenBase/EvtPatches.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

EvtItgFourCoeffFcn::EvtItgFourCoeffFcn(
    double ( *theFunction )( double, const std::vector<double>&,
                             const std::vector<double>&,
                             const std::vector<double>&,
                             const std::vector<double>& ),
    double lowerRange, double upperRange, const std::vector<double>& coeffs1,
    const std::vector<double>& coeffs2, const std::vector<double>& coeffs3,
    const std::vector<double>& coeffs4 ) :
    EvtItgAbsFunction( lowerRange, upperRange ),
    _myFunction( theFunction ),
    _coeffs1( coeffs1 ),
    _coeffs2( coeffs2 ),
    _coeffs3( coeffs3 ),
    _coeffs4( coeffs4 )
{
}

double EvtItgFourCoeffFcn::myFunction( double x ) const
{
    return _myFunction( x, _coeffs1, _coeffs2, _coeffs3, _coeffs4 );
}

void EvtItgFourCoeffFcn::setCoeff( int vect, int which, double value )
{
    if ( vect == 1 )
        _coeffs1[which] = value;
    else if ( vect == 2 )
        _coeffs2[which] = value;
    else if ( vect == 3 )
        _coeffs3[which] = value;
    else if ( vect == 4 )
        _coeffs4[which] = value;
}

double EvtItgFourCoeffFcn::getCoeff( int vect, int which )
{
    if ( vect == 1 )
        return _coeffs1[which];
    else if ( vect == 2 )
        return _coeffs2[which];
    else if ( vect == 3 )
        return _coeffs3[which];
    else if ( vect == 4 )
        return _coeffs4[which];
    else {
        return 0;
    }
}
