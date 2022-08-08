
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

#ifndef EVT_AMP_PDF_HH
#define EVT_AMP_PDF_HH

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtPdf.hh"

template <class T>

class EvtAmpPdf : public EvtPdf<T> {
  public:
    EvtAmpPdf() {}
    EvtAmpPdf( const EvtAmplitude<T>& amp ) : EvtPdf<T>(), _amp( amp.clone() )
    {
    }
    EvtAmpPdf( const EvtAmpPdf<T>& other ) :
        EvtPdf<T>( other ), COPY_PTR( _amp )
    {
    }
    virtual ~EvtAmpPdf() { delete _amp; }

    EvtAmpPdf<T>* clone() const override { return new EvtAmpPdf( *this ); }

    double pdf( const T& p ) const override
    {
        EvtComplex amp = _amp->evaluate( p );
        return real( amp ) * real( amp ) + imag( amp ) * imag( amp );
    }

  private:
    EvtAmplitude<T>* _amp;
};

#endif
