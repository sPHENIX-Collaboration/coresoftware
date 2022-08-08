
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

#ifndef EVT_VAL_ERROR_HH
#define EVT_VAL_ERROR_HH

#include <assert.h>
#include <iostream>
#include <math.h>

// Value and its associated error. E.g. this could be interval size and
// the error associated with numerical integration.

class EvtValError final {
  public:
    EvtValError();
    EvtValError( double val );
    EvtValError( double val, double err );
    EvtValError( const EvtValError& other );

    inline int valueKnown() const { return _valKnown; }
    inline double value() const
    {
        assert( _valKnown );
        return _val;
    }
    inline int errorKnown() const { return _errKnown; }
    inline double error() const
    {
        assert( _errKnown );
        return _err;
    }

    double prec() const;
    void operator=( const EvtValError& other );
    void operator*=( const EvtValError& other );
    void operator/=( const EvtValError& other );
    void operator+=( const EvtValError& other );
    void operator*=( double c );

    void print( std::ostream& ) const;

  private:
    int _valKnown;
    double _val;
    int _errKnown;
    double _err;
};

EvtValError operator*( const EvtValError& x1, const EvtValError& x2 );
EvtValError operator/( const EvtValError& x1, const EvtValError& x2 );
EvtValError operator+( const EvtValError& x1, const EvtValError& x2 );
EvtValError operator*( const EvtValError& x, double c );
EvtValError operator*( double c, const EvtValError& x );

std::ostream& operator<<( std::ostream&, const EvtValError& );

// Perform an accept/reject fraction count

template <class InputIterator, class Predicate>
EvtValError accept_reject( InputIterator it, InputIterator end, Predicate pred )
{
    int itsTried = 0;
    int itsPassed = 0;
    while ( it != end ) {
        itsTried++;
        if ( pred( *it++ ) )
            itsPassed++;
    }

    return EvtValError( ( (double)itsPassed ) / ( (double)itsTried ),
                        sqrt( itsPassed ) / ( (double)itsTried ) );
}

#endif
