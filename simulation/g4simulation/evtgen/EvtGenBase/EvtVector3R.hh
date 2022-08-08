
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

#ifndef EVTVECTOR3R_HH
#define EVTVECTOR3R_HH

#include <iosfwd>

class EvtVector3R final {
    friend EvtVector3R rotateEuler( const EvtVector3R& v, double phi,
                                    double theta, double ksi );

    inline friend EvtVector3R operator*( double c, const EvtVector3R& v2 );
    inline friend double operator*( const EvtVector3R& v1, const EvtVector3R& v2 );
    inline friend EvtVector3R operator+( const EvtVector3R& v1,
                                         const EvtVector3R& v2 );
    inline friend EvtVector3R operator-( const EvtVector3R& v1,
                                         const EvtVector3R& v2 );
    inline friend EvtVector3R operator*( const EvtVector3R& v1, double c );
    inline friend EvtVector3R operator/( const EvtVector3R& v1, double c );
    friend EvtVector3R cross( const EvtVector3R& v1, const EvtVector3R& v2 );

  public:
    EvtVector3R();
    EvtVector3R( double x, double y, double z );
    inline EvtVector3R& operator*=( const double c );
    inline EvtVector3R& operator/=( const double c );
    inline EvtVector3R& operator+=( const EvtVector3R& v2 );
    inline EvtVector3R& operator-=( const EvtVector3R& v2 );
    inline void set( int i, double d );
    inline void set( double x, double y, double z );
    void applyRotateEuler( double phi, double theta, double ksi );
    inline double get( int i ) const;
    friend std::ostream& operator<<( std::ostream& s, const EvtVector3R& v );
    double dot( const EvtVector3R& v2 );
    double d3mag() const;

  private:
    double v[3];
};

inline EvtVector3R& EvtVector3R::operator*=( const double c )
{
    v[0] *= c;
    v[1] *= c;
    v[2] *= c;
    return *this;
}

inline EvtVector3R& EvtVector3R::operator/=( const double c )
{
    v[0] /= c;
    v[1] /= c;
    v[2] /= c;
    return *this;
}

inline EvtVector3R& EvtVector3R::operator+=( const EvtVector3R& v2 )
{
    v[0] += v2.v[0];
    v[1] += v2.v[1];
    v[2] += v2.v[2];
    return *this;
}

inline EvtVector3R& EvtVector3R::operator-=( const EvtVector3R& v2 )
{
    v[0] -= v2.v[0];
    v[1] -= v2.v[1];
    v[2] -= v2.v[2];
    return *this;
}

inline EvtVector3R operator*( double c, const EvtVector3R& v2 )
{
    return EvtVector3R( v2 ) *= c;
}

inline EvtVector3R operator*( const EvtVector3R& v1, double c )
{
    return EvtVector3R( v1 ) *= c;
}

inline EvtVector3R operator/( const EvtVector3R& v1, double c )
{
    return EvtVector3R( v1 ) /= c;
}

inline double operator*( const EvtVector3R& v1, const EvtVector3R& v2 )
{
    return v1.v[0] * v2.v[0] + v1.v[1] * v2.v[1] + v1.v[2] * v2.v[2];
}

inline EvtVector3R operator+( const EvtVector3R& v1, const EvtVector3R& v2 )
{
    return EvtVector3R( v1 ) += v2;
}

inline EvtVector3R operator-( const EvtVector3R& v1, const EvtVector3R& v2 )
{
    return EvtVector3R( v1 ) -= v2;
}

inline double EvtVector3R::get( int i ) const
{
    return v[i];
}

inline void EvtVector3R::set( int i, double d )
{
    v[i] = d;
}

inline void EvtVector3R::set( double x, double y, double z )
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

#endif
