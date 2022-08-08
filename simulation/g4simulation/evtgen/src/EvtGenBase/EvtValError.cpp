
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

#include "EvtGenBase/EvtValError.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::endl;
using std::ostream;

EvtValError::EvtValError() :
    _valKnown( 0 ), _val( 0. ), _errKnown( 0 ), _err( 0. )
{
}

EvtValError::EvtValError( double val ) :
    _valKnown( 1 ), _val( val ), _errKnown( 0 ), _err( 0. )
{
}

EvtValError::EvtValError( double val, double err ) :
    _valKnown( 1 ), _val( val ), _errKnown( 1 ), _err( err )
{
}

EvtValError::EvtValError( const EvtValError& other ) :
    _valKnown( other._valKnown ),
    _val( other._val ),
    _errKnown( other._errKnown ),
    _err( other._err )
{
}

double EvtValError::prec() const
{
    assert( _valKnown && _errKnown );
    return ( _val != 0 ) ? _err / _val : 0;
}

void EvtValError::operator=( const EvtValError& other )
{
    _valKnown = other._valKnown;
    _val = other._val;
    _errKnown = other._errKnown;
    _err = other._err;
}

void EvtValError::operator*=( const EvtValError& other )
{
    assert( _valKnown && other._valKnown );

    // Relative errors add in quadrature
    if ( _errKnown && other._errKnown )
        _err = _val * other._val *
               sqrt( prec() * prec() + other.prec() * other.prec() );
    else
        _errKnown = 0;

    // Modify the value
    _val *= other._val;
}

void EvtValError::operator/=( const EvtValError& other )
{
    assert( _valKnown && other._valKnown && other._val != 0. );

    // Relative errors add in quadrature
    if ( _errKnown && other._errKnown )
        _err = _val / other._val *
               sqrt( prec() * prec() + other.prec() * other.prec() );
    else
        _errKnown = 0;

    // Modify the value
    _val /= other._val;
}

void EvtValError::print( ostream& os ) const
{
    if ( _valKnown )
        os << _val;
    else
        os << "Undef";
    os << " +/- ";
    if ( _errKnown )
        os << _err;
    else
        os << "Undef";
    os << endl;
}

void EvtValError::operator+=( const EvtValError& other )
{
    assert( _valKnown );
    assert( other._valKnown );
    _val += other._val;

    // add errors in quadrature

    if ( _errKnown && other._errKnown ) {
        _err = sqrt( _err * _err + other._err * other._err );
    } else {
        _errKnown = 0;
    }
}

void EvtValError::operator*=( double c )
{
    assert( _valKnown );
    _val *= c;
    if ( _errKnown )
        _err *= c;
}

EvtValError operator*( const EvtValError& x1, const EvtValError& x2 )
{
    EvtValError ret( x1 );
    ret *= x2;
    return ret;
}

EvtValError operator/( const EvtValError& x1, const EvtValError& x2 )
{
    EvtValError ret( x1 );
    ret /= x2;
    return ret;
}

EvtValError operator+( const EvtValError& x1, const EvtValError& x2 )
{
    EvtValError ret( x1 );
    ret += x2;
    return ret;
}

EvtValError operator*( const EvtValError& x, double c )
{
    EvtValError ret( x );
    ret *= c;
    return ret;
}

EvtValError operator*( double c, const EvtValError& x )
{
    EvtValError ret( x );
    ret *= c;
    return ret;
}

ostream& operator<<( ostream& os, const EvtValError& other )
{
    other.print( os );
    return os;
}
