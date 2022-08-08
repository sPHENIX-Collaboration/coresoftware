
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

#ifndef EVT3RANK3C_HH
#define EVT3RANK3C_HH

#include "EvtGenBase/EvtComplex.hh"

#include <iostream>

class EvtTensor3C;
class EvtVector3C;
class EvtVector3R;

class Evt3Rank3C;
inline Evt3Rank3C operator*( const EvtComplex& c, const Evt3Rank3C& t2 );
inline Evt3Rank3C operator*( const double d, const Evt3Rank3C& t2 );
inline Evt3Rank3C operator*( const Evt3Rank3C& t2, const EvtComplex& c );
inline Evt3Rank3C operator*( const Evt3Rank3C& t2, const double d );
inline Evt3Rank3C operator+( const Evt3Rank3C& t1, const Evt3Rank3C& t2 );
inline Evt3Rank3C operator-( const Evt3Rank3C& t1, const Evt3Rank3C& t2 );
Evt3Rank3C directProd( const EvtVector3C& c1, const EvtVector3C& c2,
                       const EvtVector3C& c3 );
Evt3Rank3C conj( const Evt3Rank3C& t2 );

class Evt3Rank3C final {
    friend Evt3Rank3C operator*( const EvtComplex& c, const Evt3Rank3C& t2 );
    friend Evt3Rank3C operator*( const double d, const Evt3Rank3C& t2 );
    friend Evt3Rank3C operator*( const Evt3Rank3C& t2, const EvtComplex& c );
    friend Evt3Rank3C operator*( const Evt3Rank3C& t2, const double d );
    friend Evt3Rank3C operator+( const Evt3Rank3C& t1, const Evt3Rank3C& t2 );
    friend Evt3Rank3C operator-( const Evt3Rank3C& t1, const Evt3Rank3C& t2 );
    friend Evt3Rank3C directProd( const EvtVector3C& c1, const EvtVector3C& c2,
                                  const EvtVector3C& c3 );
    friend Evt3Rank3C conj( const Evt3Rank3C& t2 );

    friend std::ostream& operator<<( std::ostream& s, const Evt3Rank3C& t2 );

  public:
    Evt3Rank3C();
    Evt3Rank3C( const Evt3Rank3C& t1 );
    Evt3Rank3C& operator=( const Evt3Rank3C& t1 );
    inline void set( int i, int j, int k, const EvtComplex& c );
    inline const EvtComplex& get( int i, int j, int k ) const;
    void zero();

    Evt3Rank3C& operator+=( const Evt3Rank3C& t2 );
    Evt3Rank3C& operator-=( const Evt3Rank3C& t2 );
    Evt3Rank3C& operator*=( const double d );
    Evt3Rank3C& operator*=( const EvtComplex& c );
    Evt3Rank3C conj() const;
    EvtTensor3C cont1( const EvtVector3C& v ) const;
    EvtTensor3C cont2( const EvtVector3C& v ) const;
    EvtTensor3C cont3( const EvtVector3C& v ) const;
    EvtTensor3C cont1( const EvtVector3R& v ) const;
    EvtTensor3C cont2( const EvtVector3R& v ) const;
    EvtTensor3C cont3( const EvtVector3R& v ) const;

  private:
    EvtComplex t[3][3][3];
};

inline Evt3Rank3C operator*( const EvtComplex& c, const Evt3Rank3C& t2 )
{
    return Evt3Rank3C( t2 ) *= c;
}

inline Evt3Rank3C operator*( const double d, const Evt3Rank3C& t2 )
{
    return Evt3Rank3C( t2 ) *= d;
}

inline Evt3Rank3C operator*( const Evt3Rank3C& t2, const EvtComplex& c )
{
    return Evt3Rank3C( t2 ) *= c;
}

inline Evt3Rank3C operator*( const Evt3Rank3C& t2, const double d )
{
    return Evt3Rank3C( t2 ) *= d;
}

inline Evt3Rank3C operator+( const Evt3Rank3C& t1, const Evt3Rank3C& t2 )
{
    return Evt3Rank3C( t1 ) += t2;
}

inline Evt3Rank3C operator-( const Evt3Rank3C& t1, const Evt3Rank3C& t2 )
{
    return Evt3Rank3C( t1 ) -= t2;
}

inline void Evt3Rank3C::set( int i, int j, int k, const EvtComplex& c )
{
    t[i][j][k] = c;
}

inline const EvtComplex& Evt3Rank3C::get( int i, int j, int k ) const
{
    return t[i][j][k];
}

#endif
