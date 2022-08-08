
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

#ifndef EVT_FLAT_AMP_HH
#define EVT_FLAT_AMP_HH

#include "EvtGenBase/EvtAmplitude.hh"

// Flat amplitude

template <class T>
class EvtFlatAmp : public EvtAmplitude<T> {
  public:
    EvtFlatAmp() {}
    EvtFlatAmp( const EvtFlatAmp<T>& other ) : EvtAmplitude<T>( other ) {}
    virtual ~EvtFlatAmp() {}

    EvtAmplitude<T>* clone() const override
    {
        return new EvtFlatAmp<T>( *this );
    }
    EvtComplex amplitude( const T& ) const override
    {
        return EvtComplex( 1., 0. );
    }
};

#endif
