
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

#ifndef EVTDIRACSPINOR_HH
#define EVTDIRACSPINOR_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4R.hh"

class EvtTensor4C;
class EvtVector4C;
class EvtDiracSpinor;

class EvtDiracSpinor final {
    friend EvtDiracSpinor rotateEuler( const EvtDiracSpinor& sp, double alpha,
                                       double beta, double gamma );
    friend EvtDiracSpinor boostTo( const EvtDiracSpinor& sp,
                                   const EvtVector4R p4 );
    friend EvtDiracSpinor boostTo( const EvtDiracSpinor& sp,
                                   const EvtVector3R boost );
    friend EvtVector4C EvtLeptonVACurrent( const EvtDiracSpinor& d,
                                           const EvtDiracSpinor& dp );
    friend EvtVector4C EvtLeptonVCurrent( const EvtDiracSpinor& d,
                                          const EvtDiracSpinor& dp );
    friend EvtVector4C EvtLeptonACurrent( const EvtDiracSpinor& d,
                                          const EvtDiracSpinor& dp );
    friend EvtComplex EvtLeptonSCurrent( const EvtDiracSpinor& d,
                                         const EvtDiracSpinor& dp );
    friend EvtComplex EvtLeptonPCurrent( const EvtDiracSpinor& d,
                                         const EvtDiracSpinor& dp );
    friend EvtTensor4C EvtLeptonTCurrent( const EvtDiracSpinor& d,
                                          const EvtDiracSpinor& dp );
    friend EvtDiracSpinor operator+( const EvtDiracSpinor& u1,
                                     const EvtDiracSpinor& u2 );
    friend EvtDiracSpinor operator-( const EvtDiracSpinor& u1,
                                     const EvtDiracSpinor& u2 );
    friend EvtDiracSpinor operator*( const EvtComplex& c,
                                     const EvtDiracSpinor& d );

    friend EvtComplex operator*( const EvtDiracSpinor& d,
                                 const EvtDiracSpinor& dp );

    friend std::ostream& operator<<( std::ostream& s, const EvtDiracSpinor& c );

  public:
    inline EvtDiracSpinor();
    EvtDiracSpinor( const EvtComplex& sp0, const EvtComplex& sp1,
                    const EvtComplex& sp2, const EvtComplex& sp3 );
    inline EvtDiracSpinor( const EvtDiracSpinor& dspinor );
    inline EvtDiracSpinor& operator=( const EvtDiracSpinor& dspinor );

    inline EvtDiracSpinor& operator+=( const EvtDiracSpinor& u2 );
    inline EvtDiracSpinor& operator-=( const EvtDiracSpinor& u2 );

    void set( const EvtComplex& sp0, const EvtComplex& sp1,
              const EvtComplex& sp2, const EvtComplex& sp3 );
    void set_spinor( int i, const EvtComplex& sp );
    const EvtComplex& get_spinor( int i ) const;
    EvtDiracSpinor conj() const;
    void applyRotateEuler( double alpha, double beta, double gamma );
    void applyBoostTo( const EvtVector4R& p4 );
    void applyBoostTo( const EvtVector3R& boost );
    EvtDiracSpinor adjoint() const;

  private:
    EvtComplex spinor[4];
};

EvtDiracSpinor::EvtDiracSpinor()
{
    spinor[0] = EvtComplex();
    spinor[1] = EvtComplex();
    spinor[2] = EvtComplex();
    spinor[3] = EvtComplex();
}

EvtDiracSpinor::EvtDiracSpinor( const EvtDiracSpinor& dspinor )
{
    spinor[0] = dspinor.spinor[0];
    spinor[1] = dspinor.spinor[1];
    spinor[2] = dspinor.spinor[2];
    spinor[3] = dspinor.spinor[3];
}

EvtDiracSpinor& EvtDiracSpinor::operator=( const EvtDiracSpinor& dspinor )
{
    spinor[0] = dspinor.spinor[0];
    spinor[1] = dspinor.spinor[1];
    spinor[2] = dspinor.spinor[2];
    spinor[3] = dspinor.spinor[3];

    return *this;
}

inline EvtDiracSpinor& EvtDiracSpinor::operator+=( const EvtDiracSpinor& u2 )
{
    spinor[0] += u2.spinor[0];
    spinor[1] += u2.spinor[1];
    spinor[2] += u2.spinor[2];
    spinor[3] += u2.spinor[3];

    return *this;
}

inline EvtDiracSpinor operator+( const EvtDiracSpinor& u1,
                                 const EvtDiracSpinor& u2 )
{
    return EvtDiracSpinor( u1 ) += u2;
}

inline EvtDiracSpinor& EvtDiracSpinor::operator-=( const EvtDiracSpinor& u2 )
{
    spinor[0] -= u2.spinor[0];
    spinor[1] -= u2.spinor[1];
    spinor[2] -= u2.spinor[2];
    spinor[3] -= u2.spinor[3];

    return *this;
}

inline EvtDiracSpinor operator-( const EvtDiracSpinor& u1,
                                 const EvtDiracSpinor& u2 )
{
    return EvtDiracSpinor( u1 ) -= u2;
}

#endif
