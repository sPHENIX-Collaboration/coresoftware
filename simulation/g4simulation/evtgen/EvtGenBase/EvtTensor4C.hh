
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

#ifndef EvtTensor4C_HH
#define EvtTensor4C_HH

#include "EvtGenBase/EvtComplex.hh"

//Class to handle 4D complex valued tensors.
class EvtTensor4C;
class EvtVector4C;
class EvtVector4R;
class EvtVector3R;

namespace EvtGenFunctions {
    EvtTensor4C directProd( const EvtVector4R& c1, const EvtVector4R& c2 );
    EvtTensor4C directProd( const EvtVector4C& c1, const EvtVector4C& c2 );
    EvtTensor4C directProd( const EvtVector4C& c1, const EvtVector4R& c2 );
}    // namespace EvtGenFunctions

class EvtTensor4C final {
    friend EvtTensor4C EvtGenFunctions::directProd( const EvtVector4R& c1,
                                                    const EvtVector4R& c2 );
    friend EvtTensor4C EvtGenFunctions::directProd( const EvtVector4C& c1,
                                                    const EvtVector4C& c2 );
    friend EvtTensor4C EvtGenFunctions::directProd( const EvtVector4C& c1,
                                                    const EvtVector4R& c2 );

    friend EvtTensor4C rotateEuler( const EvtTensor4C& e, double alpha,
                                    double beta, double gamma );
    friend EvtTensor4C boostTo( const EvtTensor4C& e, const EvtVector4R p4 );
    friend EvtTensor4C boostTo( const EvtTensor4C& e, const EvtVector3R boost );
    friend EvtTensor4C dual( const EvtTensor4C& t2 );
    friend EvtTensor4C conj( const EvtTensor4C& t2 );
    friend EvtTensor4C cont22( const EvtTensor4C& t1, const EvtTensor4C& t2 );
    friend EvtTensor4C cont11( const EvtTensor4C& t1, const EvtTensor4C& t2 );
    friend EvtTensor4C operator*( const EvtTensor4C& t1, const EvtComplex& c );
    friend EvtTensor4C operator*( const EvtComplex& c, const EvtTensor4C& t1 );
    friend EvtTensor4C operator*( const EvtTensor4C& t1, double d );
    friend EvtTensor4C operator*( double d, const EvtTensor4C& t1 );
    friend EvtComplex cont( const EvtTensor4C& t1, const EvtTensor4C& t2 );
    friend EvtTensor4C operator+( const EvtTensor4C& t1, const EvtTensor4C& t2 );
    friend EvtTensor4C operator-( const EvtTensor4C& t1, const EvtTensor4C& t2 );

  public:
    EvtTensor4C() { ; }

    EvtTensor4C( double t00, double t11, double t22, double t33 )
    {
        setdiag( t00, t11, t22, t33 );
    }

    EvtTensor4C( const EvtTensor4C& t1 );
    EvtTensor4C& operator=( const EvtTensor4C& t1 );
    EvtTensor4C& operator*=( const EvtComplex& c );
    EvtTensor4C& operator*=( double d );
    EvtTensor4C& addDirProd( const EvtVector4R& p1, const EvtVector4R& p2 );
    static const EvtTensor4C& g();
    inline void set( int i, int j, const EvtComplex& c );
    void setdiag( double t00, double t11, double t22, double t33 );
    inline const EvtComplex& get( int i, int j ) const;
    inline EvtComplex trace() const;
    void zero();
    void applyRotateEuler( double alpha, double beta, double gamma );
    void applyBoostTo( const EvtVector4R& p4 );
    void applyBoostTo( const EvtVector3R& boost );
    friend std::ostream& operator<<( std::ostream& s, const EvtTensor4C& t );
    EvtTensor4C& operator+=( const EvtTensor4C& t2 );
    EvtTensor4C& operator-=( const EvtTensor4C& t2 );
    EvtTensor4C conj() const;
    EvtVector4C cont1( const EvtVector4C& v4 ) const;
    EvtVector4C cont2( const EvtVector4C& v4 ) const;
    EvtVector4C cont1( const EvtVector4R& v4 ) const;
    EvtVector4C cont2( const EvtVector4R& v4 ) const;

  private:
    EvtComplex t[4][4];
};

inline EvtTensor4C operator+( const EvtTensor4C& t1, const EvtTensor4C& t2 )
{
    return EvtTensor4C( t1 ) += t2;
}

inline EvtTensor4C operator-( const EvtTensor4C& t1, const EvtTensor4C& t2 )
{
    return EvtTensor4C( t1 ) -= t2;
}

inline void EvtTensor4C::set( int i, int j, const EvtComplex& c )
{
    t[i][j] = c;
}

inline const EvtComplex& EvtTensor4C::get( int i, int j ) const
{
    return t[i][j];
}

inline EvtComplex EvtTensor4C::trace() const
{
    return t[0][0] - t[1][1] - t[2][2] - t[3][3];
}

#endif
