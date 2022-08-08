
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

#include "EvtGenBase/EvtPropFlatte.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <iostream>
#include <math.h>
using std::cout;
using std::endl;

EvtPropFlatte::EvtPropFlatte( double m0, double g0, double m0a, double m0b,
                              double g1, double m1a, double m1b ) :
    EvtPropagator( m0, g0 ),
    _m0a( m0a ),
    _m0b( m0b ),
    _g1( g1 ),
    _m1a( m1a ),
    _m1b( m1b )
{
}

EvtAmplitude<EvtPoint1D>* EvtPropFlatte::clone() const
{
    return new EvtPropFlatte( *this );
}

EvtComplex EvtPropFlatte::amplitude( const EvtPoint1D& x ) const
{
    /*

  Use BES parameterization:

                        1.
      -----------------------------------------
       m0^2 - m^2 - i*m0*( g1*rho1 + g2*rho2 )


  Resonance mass: m0
  Channel1: m0a, m0b, g0
  Channel2: m1a, m1b, g1

  where breakup momenta q's are:

          E0a = (m^2 + m0a^2 - m0b^2) / 2m
          q0  = sqrt( E0a^2 - m0a^2 )

          E1a = (m^2 + m1a^2 - m1b^2) / 2m
          q1  = sqrt( E1a^2 - m1a^2 )


  */

    double s = x.value() * x.value();
    double m = x.value();

    double E0a = 0.5 * ( s + _m0a * _m0a - _m0b * _m0b ) / m;
    double qSq0 = E0a * E0a - _m0a * _m0a;

    double E1a = 0.5 * ( s + _m1a * _m1a - _m1b * _m1b ) / m;
    double qSq1 = E1a * E1a - _m1a * _m1a;

    EvtComplex gamma0 = qSq0 >= 0 ? EvtComplex( _g0 * sqrt( qSq0 ), 0 )
                                  : EvtComplex( 0, _g0 * sqrt( -qSq0 ) );
    EvtComplex gamma1 = qSq1 >= 0 ? EvtComplex( _g1 * sqrt( qSq1 ), 0 )
                                  : EvtComplex( 0, _g1 * sqrt( -qSq1 ) );

    EvtComplex gamma = gamma0 + gamma1;

    EvtComplex a = 1.0 /
                   ( _m0 * _m0 - s - EvtComplex( 0.0, 2 * _m0 / m ) * gamma );

    return a;
}
