
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

#ifndef EVTFLATTE_HH
#define EVTFLATTE_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <vector>

using std::vector;

// Helper class

class EvtFlatteParam {
  public:
    EvtFlatteParam( double m1, double m2, double g ) :
        _m1( m1 ), _m2( m2 ), _g( g )
    {
    }

    inline double m1() const { return _m1; }
    inline double m2() const { return _m2; }
    inline double g() const { return _g; }

  private:
    double _m1, _m2, _g;
};

//class declaration

class EvtFlatte final {
  public:
    //operator
    EvtFlatte& operator=( const EvtFlatte& );

    //constructor with all information about the resonance
    EvtFlatte( const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
               const EvtVector4R& p4_d2, double ampl, double theta, double mass,
               vector<EvtFlatteParam>& params
               //           double m1a = 0.0, double m1b = 0.0, double g1 = 0.0,
               //           double m2a = 0.0, double m2b = 0.0, double g2 = 0.0
    );

    //accessors
    //return 4-momenta of the particles involved
    inline const EvtVector4R& p4_p() { return _p4_p; }
    inline const EvtVector4R& p4_d1() { return _p4_d1; }
    inline const EvtVector4R& p4_d2() { return _p4_d2; }

    //return amplitude
    inline double amplitude() { return _ampl; }

    //return theta
    inline double theta() { return _theta; }

    //return bwm
    inline double mass() { return _mass; }

    //functions

    //calculate amplitude for this resonance
    EvtComplex resAmpl();

  private:
    inline EvtComplex sqrtCplx( double in )
    {
        return ( in > 0 ) ? EvtComplex( sqrt( in ), 0 )
                          : EvtComplex( 0, sqrt( -in ) );
    }

    EvtVector4R _p4_p, _p4_d1, _p4_d2;
    double _ampl, _theta, _mass;
    vector<EvtFlatteParam> _params;
    //      double _m1a, _m1b, _g1;
    //      double _m2a, _m2b, _g2;
};

#endif
