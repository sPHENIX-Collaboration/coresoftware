
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

#ifndef EVTRESONANCE2_HH
#define EVTRESONANCE2_HH

#include "EvtGenBase/EvtVector4R.hh"

class EvtComplex;

class EvtResonance2 {
  public:
    //operator
    EvtResonance2& operator=( const EvtResonance2& );

    // Constructor
    EvtResonance2( const EvtVector4R& p4_p, const EvtVector4R& p4_d1,
                   const EvtVector4R& p4_d2, double ampl = 1.0,
                   double theta = 0.0, double gamma = 0.0, double bwm = 0.0,
                   int spin = 0, bool invmass_angdenom = false,
                   double barrier1 = 1.5, double barrier2 = 5.0 );

    //accessors
    //return 4-momenta of the particles involved
    inline const EvtVector4R& p4_p() const { return _p4_p; }
    inline const EvtVector4R& p4_d1() const { return _p4_d1; }
    inline const EvtVector4R& p4_d2() const { return _p4_d2; }

    //return amplitude
    inline double amplitude() const { return _ampl; }

    //return theta
    inline double theta() const { return _theta; }

    //return gamma
    inline double gamma() const { return _gamma; }

    //return bwm
    inline double bwm() const { return _bwm; }

    //return spin
    inline int spin() const { return _spin; }

    //calculate amplitude for this resonance
    EvtComplex resAmpl() const;

  private:
    EvtVector4R _p4_p, _p4_d1, _p4_d2;
    double _ampl, _theta, _gamma, _bwm, _barrier1, _barrier2;
    int _spin;
    bool _invmass_angdenom;
};

#endif
