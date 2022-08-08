
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

#include "EvtGenBase/EvtComplex.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <iostream>
#include <math.h>
using std::ostream;

ostream& operator<<( ostream& s, const EvtComplex& c )
{
    s << "(" << c._rpart << "," << c._ipart << ")";
    return s;
}

EvtComplex& EvtComplex::operator*=( EvtComplex c )
{
    double r = _rpart * c._rpart - _ipart * c._ipart;
    double i = _rpart * c._ipart + _ipart * c._rpart;

    _rpart = r;
    _ipart = i;

    return *this;
}

EvtComplex& EvtComplex::operator/=( EvtComplex c )
{
    double inv = 1.0 / ( c._rpart * c._rpart + c._ipart * c._ipart );

    double r = inv * ( _rpart * c._rpart + _ipart * c._ipart );
    double i = inv * ( _ipart * c._rpart - _rpart * c._ipart );

    _rpart = r;
    _ipart = i;

    return *this;
}
