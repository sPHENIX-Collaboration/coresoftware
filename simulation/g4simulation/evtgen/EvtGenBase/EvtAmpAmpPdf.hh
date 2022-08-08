
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

#ifndef EVT_AMP_AMP_PDF_HH
#define EVT_AMP_AMP_PDF_HH

// From the product A1A2* four PDF terms can be constructed, by taking the positive
// and the negative parts or the real and imaginary part of the product.

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtPdf.hh"

#include <assert.h>

enum
{
    POSRE = 0,
    NEGRE,
    POSIM,
    NEGIM
};

template <class T>
class EvtAmpAmpPdf : public EvtPdf<T> {
  public:
    EvtAmpAmpPdf() {}
    EvtAmpAmpPdf( int type, const EvtAmplitude<T>& amp1,
                  const EvtAmplitude<T>& amp2 ) :
        EvtPdf<T>(), _type( type ), _amp1( amp1.clone() ), _amp2( amp2.clone() )
    {
    }
    EvtAmpAmpPdf( const EvtAmpAmpPdf<T>& other ) :
        EvtPdf<T>( other ),
        _type( other._type ),
        COPY_PTR( _amp1 ),
        COPY_PTR( _amp2 )
    {
    }
    virtual ~EvtAmpAmpPdf()
    {
        delete _amp1;
        delete _amp2;
    }

    virtual EvtAmpAmpPdf<T>* clone() const { return new EvtAmpAmpPdf( *this ); }

    virtual double pdf( const T& p ) const
    {
        EvtComplex amp1 = _amp1->evaluate( p );
        EvtComplex amp2 = _amp2->evaluate( p );
        EvtComplex pr = amp1 * conj( amp2 );

        if ( _type == POSRE )
            return real( pr ) > 0 ? real( pr ) : 0.;
        if ( _type == NEGRE )
            return real( pr ) < 0 ? -real( pr ) : 0.;
        if ( _type == POSIM )
            return imag( pr ) > 0 ? imag( pr ) : 0.;
        if ( _type == NEGIM )
            return imag( pr ) < 0 ? -imag( pr ) : 0.;

        assert( 0 );
    }

  private:
    int _type;
    EvtAmplitude<T>* _amp1;
    EvtAmplitude<T>* _amp2;
};

#endif
