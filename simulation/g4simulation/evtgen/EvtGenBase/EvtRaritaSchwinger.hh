
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

#ifndef EVTRARITASCHWINGER_HH
#define EVTRARITASCHWINGER_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

class EvtRaritaSchwinger;
EvtRaritaSchwinger rotateEuler( const EvtRaritaSchwinger& rs, double alpha,
                                double beta, double gamma );
EvtRaritaSchwinger boostTo( const EvtRaritaSchwinger& rs, const EvtVector4R p4 );
EvtRaritaSchwinger boostTo( const EvtRaritaSchwinger& rs,
                            const EvtVector3R boost );
EvtRaritaSchwinger dirProd( EvtVector4R v, EvtDiracSpinor u );
EvtRaritaSchwinger dirProd( EvtVector4C v, EvtDiracSpinor u );
EvtRaritaSchwinger operator+( const EvtRaritaSchwinger& u1,
                              const EvtRaritaSchwinger& u2 );
EvtRaritaSchwinger operator-( const EvtRaritaSchwinger& u1,
                              const EvtRaritaSchwinger& u2 );
EvtComplex operator*( const EvtRaritaSchwinger& u1, const EvtRaritaSchwinger& u2 );

class EvtRaritaSchwinger final {
    friend EvtRaritaSchwinger rotateEuler( const EvtRaritaSchwinger& rs,
                                           double alpha, double beta,
                                           double gamma );
    friend EvtRaritaSchwinger boostTo( const EvtRaritaSchwinger& rs,
                                       const EvtVector4R p4 );
    friend EvtRaritaSchwinger boostTo( const EvtRaritaSchwinger& rs,
                                       const EvtVector3R boost );

    friend EvtRaritaSchwinger dirProd( EvtVector4R v, EvtDiracSpinor u );
    friend EvtRaritaSchwinger dirProd( EvtVector4C v, EvtDiracSpinor u );

    friend EvtRaritaSchwinger operator+( const EvtRaritaSchwinger& u1,
                                         const EvtRaritaSchwinger& u2 );
    friend EvtRaritaSchwinger operator-( const EvtRaritaSchwinger& u1,
                                         const EvtRaritaSchwinger& u2 );

    friend EvtComplex operator*( const EvtRaritaSchwinger& u1,
                                 const EvtRaritaSchwinger& u2 );

  public:
    inline EvtRaritaSchwinger();
    inline EvtRaritaSchwinger( const EvtRaritaSchwinger& rs );
    inline EvtRaritaSchwinger& operator=( const EvtRaritaSchwinger& rs );

    void set( int i, int j, const EvtComplex& sp );

    void applyRotateEuler( double alpha, double beta, double gamma );
    void applyBoostTo( const EvtVector4R p4 );
    void applyBoostTo( const EvtVector3R boost );

    EvtRaritaSchwinger& operator+=( const EvtRaritaSchwinger& u2 );
    EvtRaritaSchwinger& operator-=( const EvtRaritaSchwinger& u2 );

    EvtComplex get( int i, int j ) const;
    friend std::ostream& operator<<( std::ostream& s,
                                     const EvtRaritaSchwinger& rs );

    EvtVector4C getVector( int i ) const;
    EvtDiracSpinor getSpinor( int i ) const;

    void setVector( int i, const EvtVector4C& v );
    void setSpinor( int i, const EvtDiracSpinor& sp );

  private:
    //First index in spinor index, second is Lorentz index.
    EvtComplex _rs[4][4];
};

EvtRaritaSchwinger::EvtRaritaSchwinger()
{
    int i, j;
    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            _rs[i][j] = 0.0;
        }
    }
}

EvtRaritaSchwinger::EvtRaritaSchwinger( const EvtRaritaSchwinger& rs )
{
    int i, j;
    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            _rs[i][j] = rs._rs[i][j];
        }
    }
}

EvtRaritaSchwinger& EvtRaritaSchwinger::operator=( const EvtRaritaSchwinger& rs )
{
    int i, j;
    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ) {
            _rs[i][j] = rs._rs[i][j];
        }
    }

    return *this;
}

#endif
