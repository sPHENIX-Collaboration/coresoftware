
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

#ifndef EVTGAMMAMATRIX_HH
#define EVTGAMMAMATRIX_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"    // needed for adjoint
//#include <iostream.h>
#include <iosfwd>
class EvtGammaMatrix;
class EvtVector4C;

namespace EvtGenFunctions {
    // slash or Feynman slash a 4-vector
    EvtGammaMatrix slash( const EvtVector4C& p );
    EvtGammaMatrix slash( const EvtVector4R& p );
}    // namespace EvtGenFunctions

// Description: Class to manipulate gamma matrices. The reperesentation
//              used is the "standard" Dirac representation.

class EvtGammaMatrix final {
    friend EvtGammaMatrix operator*( const EvtComplex& c,
                                     const EvtGammaMatrix& g );
    friend EvtGammaMatrix operator*( const EvtGammaMatrix& g,
                                     const EvtComplex& c );
    friend EvtGammaMatrix operator/( const EvtGammaMatrix& g, const double d );
    friend EvtDiracSpinor operator*( const EvtGammaMatrix& g,
                                     const EvtDiracSpinor& d );
    friend EvtGammaMatrix operator+( const EvtGammaMatrix& g1,
                                     const EvtGammaMatrix& g2 );
    friend EvtGammaMatrix operator-( const EvtGammaMatrix& g1,
                                     const EvtGammaMatrix& g2 );
    friend EvtGammaMatrix operator*( const EvtGammaMatrix& g1,
                                     const EvtGammaMatrix& g2 );
    friend std::ostream& operator<<( std::ostream& s, const EvtGammaMatrix& v );
    friend EvtDiracSpinor EvtDiracSpinor::adjoint() const;

  public:
    EvtGammaMatrix();
    EvtGammaMatrix( const EvtGammaMatrix& gm );
    EvtGammaMatrix& operator=( const EvtGammaMatrix& gm );

    void init();
    static const EvtGammaMatrix& g( int );
    static const EvtGammaMatrix& g0();
    static const EvtGammaMatrix& g1();
    static const EvtGammaMatrix& g2();
    static const EvtGammaMatrix& g3();
    static const EvtGammaMatrix& g5();
    static const EvtGammaMatrix& id();
    static const EvtGammaMatrix& va0();
    static const EvtGammaMatrix& va1();
    static const EvtGammaMatrix& va2();
    static const EvtGammaMatrix& va3();
    static const EvtGammaMatrix& v0();
    static const EvtGammaMatrix& v1();
    static const EvtGammaMatrix& v2();
    static const EvtGammaMatrix& v3();
    // Dirac sigma matrix with upper or lower indices (only one element)
    static const EvtGammaMatrix& sigmaUpper( unsigned int mu, unsigned int nu );
    static const EvtGammaMatrix& sigmaLower( unsigned int mu, unsigned int nu );

    EvtGammaMatrix& operator+=( const EvtGammaMatrix& g );
    EvtGammaMatrix& operator-=( const EvtGammaMatrix& g );
    EvtGammaMatrix& operator*=( const EvtGammaMatrix& g );

  private:
    EvtComplex _gamma[4][4];
};

inline EvtGammaMatrix operator+( const EvtGammaMatrix& g1,
                                 const EvtGammaMatrix& g2 )
{
    return EvtGammaMatrix( g1 ) += g2;
}

inline EvtGammaMatrix operator-( const EvtGammaMatrix& g1,
                                 const EvtGammaMatrix& g2 )
{
    return EvtGammaMatrix( g1 ) -= g2;
}

inline EvtGammaMatrix operator*( const EvtGammaMatrix& g1,
                                 const EvtGammaMatrix& g2 )
{
    return EvtGammaMatrix( g1 ) *= g2;
}

inline EvtGammaMatrix operator/( const EvtGammaMatrix& g, const double d )
{
    return g * EvtComplex( 1 / d, 0 );
}

#endif
