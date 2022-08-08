
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

#ifndef EVT_PROPAGATOR_HH
#define EVT_PROPAGATOR_HH

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPoint1D.hh"

#include <assert.h>

// Defines propagator as a function of mass and width

class EvtPropagator : public EvtAmplitude<EvtPoint1D> {
  public:
    EvtPropagator( double m0, double g0 ) : _m0( m0 ), _g0( g0 )
    {
        assert( m0 > 0 );
        assert( g0 >= 0 );
    }

    // Accessors

    inline double m0() const { return _m0; }
    inline double g0() const { return _g0; }

    // Modifiers (can be useful e.g. for fitting!)

    inline void set_m0( double m0 )
    {
        assert( m0 > 0 );
        _m0 = m0;
    }
    inline void set_g0( double g0 )
    {
        assert( g0 >= 0 );
        _g0 = g0;
    }

  protected:
    double _m0;
    double _g0;
};

#endif
