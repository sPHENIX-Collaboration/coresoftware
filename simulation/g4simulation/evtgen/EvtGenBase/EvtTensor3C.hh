
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

#ifndef EVTTENSOR3C_HH
#define EVTTENSOR3C_HH

#include "EvtGenBase/EvtComplex.hh"

#include <iostream>

class EvtVector3C;
class EvtVector3R;

class EvtTensor3C;

namespace EvtGenFunctions {
    EvtTensor3C eps( const EvtVector3R& v );
    EvtTensor3C rotateEuler( const EvtTensor3C& v, double phi, double theta,
                             double ksi );
    EvtTensor3C directProd( const EvtVector3C& c1, const EvtVector3C& c2 );
    EvtTensor3C directProd( const EvtVector3C& c1, const EvtVector3R& c2 );
    EvtTensor3C directProd( const EvtVector3R& c1, const EvtVector3R& c2 );
}    // namespace EvtGenFunctions

class EvtTensor3C final {
    friend EvtTensor3C operator*( const EvtComplex& c, const EvtTensor3C& t2 );
    friend EvtTensor3C operator*( const double d, const EvtTensor3C& t2 );
    friend EvtTensor3C operator*( const EvtTensor3C& t2, const EvtComplex& c );
    friend EvtTensor3C operator*( const EvtTensor3C& t2, const double d );
    friend EvtTensor3C operator+( const EvtTensor3C& t1, const EvtTensor3C& t2 );
    friend EvtTensor3C operator-( const EvtTensor3C& t1, const EvtTensor3C& t2 );
    friend EvtTensor3C EvtGenFunctions::directProd( const EvtVector3C& c1,
                                                    const EvtVector3C& c2 );
    friend EvtTensor3C EvtGenFunctions::directProd( const EvtVector3C& c1,
                                                    const EvtVector3R& c2 );
    friend EvtTensor3C EvtGenFunctions::directProd( const EvtVector3R& c1,
                                                    const EvtVector3R& c2 );
    friend EvtTensor3C conj( const EvtTensor3C& t2 );
    //Contract the second index of two tensors result(i,j) = t1(i,k)t2(j,k)
    friend EvtTensor3C cont22( const EvtTensor3C& t1, const EvtTensor3C& t2 );
    //Contract the first index of two tensors result(i,j) = t1(k,i)t2(k,j)
    friend EvtTensor3C cont11( const EvtTensor3C& t1, const EvtTensor3C& t2 );
    //Contract the last index of eps_{ijk} with w
    friend EvtTensor3C EvtGenFunctions::eps( const EvtVector3R& v );
    friend std::ostream& operator<<( std::ostream& c, const EvtTensor3C& v );

  public:
    EvtTensor3C();
    EvtTensor3C( const EvtTensor3C& t1 );
    EvtTensor3C( double d11, double d22, double d33 );
    EvtTensor3C& operator=( const EvtTensor3C& t1 );
    inline void set( int i, int j, const EvtComplex& c );
    inline const EvtComplex& get( int i, int j ) const;
    inline EvtComplex trace() const;
    static const EvtTensor3C& id();
    void zero();
    void applyRotateEuler( double phi, double theta, double ksi );

    EvtTensor3C operator+=( const EvtTensor3C& t2 );
    EvtTensor3C operator-=( const EvtTensor3C& t2 );
    EvtTensor3C operator*=( const double d );
    EvtTensor3C operator*=( const EvtComplex& c );
    EvtTensor3C conj() const;
    EvtVector3C cont1( const EvtVector3C& v ) const;
    EvtVector3C cont2( const EvtVector3C& v ) const;
    EvtVector3C cont1( const EvtVector3R& v ) const;
    EvtVector3C cont2( const EvtVector3R& v ) const;

  private:
    EvtComplex t[3][3];
};

inline EvtTensor3C operator*( const EvtComplex& c, const EvtTensor3C& t2 )
{
    return EvtTensor3C( t2 ) *= c;
}

inline EvtTensor3C operator*( const double d, const EvtTensor3C& t2 )
{
    return EvtTensor3C( t2 ) *= d;
}

inline EvtTensor3C operator*( const EvtTensor3C& t2, const EvtComplex& c )
{
    return EvtTensor3C( t2 ) *= c;
}

inline EvtTensor3C operator*( const EvtTensor3C& t2, const double d )
{
    return EvtTensor3C( t2 ) *= d;
}

inline EvtTensor3C operator+( const EvtTensor3C& t1, const EvtTensor3C& t2 )
{
    return EvtTensor3C( t1 ) += t2;
}

inline EvtTensor3C operator-( const EvtTensor3C& t1, const EvtTensor3C& t2 )
{
    return EvtTensor3C( t1 ) -= t2;
}

inline void EvtTensor3C::set( int i, int j, const EvtComplex& c )
{
    t[i][j] = c;
}

inline const EvtComplex& EvtTensor3C::get( int i, int j ) const
{
    return t[i][j];
}

inline EvtComplex EvtTensor3C::trace() const
{
    return t[0][0] + t[1][1] + t[2][2];
}

#endif
