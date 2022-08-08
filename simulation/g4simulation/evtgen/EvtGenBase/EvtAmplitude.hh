
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

#ifndef EVT_AMPLITUDE_HH
#define EVT_AMPLITUDE_HH

#include "EvtGenBase/EvtComplex.hh"

// Complex-valued amplitude

template <class T>
class EvtAmplitude {
  public:
    EvtAmplitude() = default;
    EvtAmplitude( const EvtAmplitude& ) = default;
    EvtAmplitude( EvtAmplitude&& ) = default;
    EvtAmplitude& operator=( const EvtAmplitude& ) = default;
    EvtAmplitude& operator=( EvtAmplitude&& ) = default;
    virtual ~EvtAmplitude() = default;

    virtual EvtAmplitude<T>* clone() const = 0;

    EvtComplex evaluate( const T& p ) const
    {
        EvtComplex ret( 0., 0. );
        if ( p.isValid() )
            ret = amplitude( p );
        return ret;
    }

  protected:
    // Derive in subclasses to define amplitude computation
    // for a fully constructed amplitude object.

    virtual EvtComplex amplitude( const T& ) const = 0;
};

#endif
