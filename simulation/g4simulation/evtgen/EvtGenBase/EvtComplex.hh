
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

#ifndef EVTCOMPLEX_HH
#define EVTCOMPLEX_HH

#include "EvtGenBase/EvtConst.hh"

#include <iostream>
#include <math.h>

class EvtComplex {
    inline friend EvtComplex operator*( double d, const EvtComplex& c );
    inline friend EvtComplex operator*( const EvtComplex& c, double d );
    inline friend EvtComplex operator/( const EvtComplex& c, double d );
    inline friend EvtComplex operator/( double d, const EvtComplex& c );
    inline friend EvtComplex operator*( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator/( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator+( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator-( const EvtComplex& c1,
                                        const EvtComplex& c2 );
    inline friend EvtComplex operator-( const EvtComplex& c );
    inline friend EvtComplex conj( const EvtComplex& c );
    inline friend double abs( const EvtComplex& c );
    inline friend double abs2( const EvtComplex& c );
    inline friend double arg( const EvtComplex& c );
    inline friend double real( const EvtComplex& c );
    inline friend double imag( const EvtComplex& c );
    inline friend EvtComplex exp( const EvtComplex& c );
    friend std::ostream& operator<<( std::ostream& s, const EvtComplex& c );

  public:
    EvtComplex() : _rpart( 0.0 ), _ipart( 0.0 ) {}
    EvtComplex( double rpart, double ipart = 0.0 ) :
        _rpart( rpart ), _ipart( ipart )
    {
    }
    EvtComplex( const EvtComplex& c ) : _rpart( c._rpart ), _ipart( c._ipart )
    {
    }
    inline EvtComplex& operator*=( double d );
    inline EvtComplex& operator/=( double d );
    EvtComplex& operator*=( EvtComplex c );
    EvtComplex& operator/=( EvtComplex c );
    inline EvtComplex& operator=( const EvtComplex& c );
    inline EvtComplex& operator+=( const EvtComplex& c );
    inline EvtComplex& operator-=( const EvtComplex& c );
    inline EvtComplex& operator+=( double d );
    inline EvtComplex& operator-=( double d );
    inline int operator==( const EvtComplex c );
    inline int operator!=( const EvtComplex c );

  private:
    double _rpart, _ipart;
};

typedef EvtComplex* EvtComplexPtr;
typedef EvtComplexPtr* EvtComplexPtrPtr;
typedef EvtComplexPtrPtr* EvtComplexPtrPtrPtr;

EvtComplex& EvtComplex::operator=( const EvtComplex& c )
{
    _rpart = c._rpart;
    _ipart = c._ipart;

    return *this;
}

EvtComplex& EvtComplex::operator+=( const EvtComplex& c )
{
    _rpart += c._rpart;
    _ipart += c._ipart;

    return *this;
}

EvtComplex& EvtComplex::operator-=( const EvtComplex& c )
{
    _rpart -= c._rpart;
    _ipart -= c._ipart;

    return *this;
}

EvtComplex& EvtComplex::operator+=( double d )
{
    _rpart += d;

    return *this;
}

EvtComplex& EvtComplex::operator-=( double d )
{
    _rpart -= d;

    return *this;
}

EvtComplex operator*( double d, const EvtComplex& c )
{
    return EvtComplex( c._rpart * d, c._ipart * d );
}

EvtComplex operator*( const EvtComplex& c, double d )
{
    return EvtComplex( c._rpart * d, c._ipart * d );
}

EvtComplex operator/( const EvtComplex& c, double d )
{
    return EvtComplex( c._rpart / d, c._ipart / d );
}

EvtComplex& EvtComplex::operator*=( double d )
{
    _rpart *= d;
    _ipart *= d;

    return *this;
}

EvtComplex& EvtComplex::operator/=( double d )
{
    _rpart /= d;
    _ipart /= d;

    return *this;
}

EvtComplex operator/( double d, const EvtComplex& c )
{
    double Num = d / ( c._rpart * c._rpart + c._ipart * c._ipart );

    return EvtComplex( Num * c._rpart, -Num * c._ipart );
}

EvtComplex operator/( const EvtComplex& c1, const EvtComplex& c2 )
{
    double inv = 1.0 / ( c2._rpart * c2._rpart + c2._ipart * c2._ipart );

    return EvtComplex( inv * ( c1._rpart * c2._rpart + c1._ipart * c2._ipart ),
                       inv * ( c1._ipart * c2._rpart - c1._rpart * c2._ipart ) );
}

EvtComplex operator*( const EvtComplex& c1, const EvtComplex& c2 )
{
    return EvtComplex( c1._rpart * c2._rpart - c1._ipart * c2._ipart,
                       c1._rpart * c2._ipart + c1._ipart * c2._rpart );
}

EvtComplex operator-( const EvtComplex& c1, const EvtComplex& c2 )
{
    return EvtComplex( c1._rpart - c2._rpart, c1._ipart - c2._ipart );
}

EvtComplex operator+( const EvtComplex& c1, const EvtComplex& c2 )
{
    return EvtComplex( c1._rpart + c2._rpart, c1._ipart + c2._ipart );
}

int EvtComplex::operator==( const EvtComplex c )
{
    return _rpart == c._rpart && _ipart == c._ipart;
}

int EvtComplex::operator!=( const EvtComplex c )
{
    return _rpart != c._rpart || _ipart != c._ipart;
}

EvtComplex operator-( const EvtComplex& c )
{
    return EvtComplex( -c._rpart, -c._ipart );
}

EvtComplex conj( const EvtComplex& c )
{
    return EvtComplex( c._rpart, -c._ipart );
}

double abs( const EvtComplex& c )
{
    double c2 = c._rpart * c._rpart + c._ipart * c._ipart;
    if ( c2 <= 0.0 )
        return 0.0;
    return sqrt( c2 );
}

double abs2( const EvtComplex& c )
{
    return c._rpart * c._rpart + c._ipart * c._ipart;
}

double arg( const EvtComplex& c )
{
    if ( ( c._rpart == 0 ) && ( c._ipart == 0 ) ) {
        return 0.0;
    }
    if ( c._rpart == 0 ) {
        if ( c._ipart > 0 ) {
            return EvtConst::pi / 2;
        } else {
            return -EvtConst::pi / 2;
        }
    } else {
        return atan2( c._ipart, c._rpart );
    }
}

double real( const EvtComplex& c )
{
    return c._rpart;
}

double imag( const EvtComplex& c )
{
    return c._ipart;
}

EvtComplex exp( const EvtComplex& c )
{
    return exp( c._rpart ) * EvtComplex( cos( c._ipart ), sin( c._ipart ) );
}

#endif
