
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

#ifndef EVT_PDF_SUM_HH
#define EVT_PDF_SUM_HH

#include <stdio.h>
#include <vector>
using std::vector;
#include "EvtGenBase/EvtPdf.hh"

// Sum of PDF functions.

template <class T>
class EvtPdfSum : public EvtPdf<T> {
  public:
    EvtPdfSum() {}
    EvtPdfSum( const EvtPdfSum<T>& other );
    virtual ~EvtPdfSum();
    EvtPdfSum* clone() const override { return new EvtPdfSum( *this ); }

    // Manipulate terms and coefficients

    void addTerm( double c, const EvtPdf<T>& pdf )
    {
        assert( c >= 0. );
        _c.push_back( c );
        _term.push_back( pdf.clone() );
    }

    void addOwnedTerm( double c, std::unique_ptr<EvtPdf<T>> pdf )
    {
        _c.push_back( c );
        _term.push_back( pdf.release() );
    }

    size_t nTerms() const { return _term.size(); }    // number of terms

    inline double c( int i ) const { return _c[i]; }
    inline EvtPdf<T>* getPdf( int i ) const { return _term[i]; }

    // Integrals

    EvtValError compute_integral() const override;
    EvtValError compute_integral( int N ) const override;
    T randomPoint() override;

  protected:
    double pdf( const T& p ) const override;

    vector<double> _c;           // coefficients
    vector<EvtPdf<T>*> _term;    // pointers to pdfs
};

template <class T>
EvtPdfSum<T>::EvtPdfSum( const EvtPdfSum<T>& other ) : EvtPdf<T>( other )
{
    for ( size_t i = 0; i < other.nTerms(); i++ ) {
        _c.push_back( other._c[i] );
        _term.push_back( other._term[i]->clone() );
    }
}

template <class T>
EvtPdfSum<T>::~EvtPdfSum()
{
    for ( size_t i = 0; i < _c.size(); i++ ) {
        delete _term[i];
    }
}

template <class T>
double EvtPdfSum<T>::pdf( const T& p ) const
{
    double ret = 0.;
    for ( size_t i = 0; i < _c.size(); i++ ) {
        ret += _c[i] * _term[i]->evaluate( p );
    }
    return ret;
}

/*
 * Compute the sum integral by summing all term integrals.
 */

template <class T>
EvtValError EvtPdfSum<T>::compute_integral() const
{
    EvtValError itg( 0.0, 0.0 );
    for ( size_t i = 0; i < nTerms(); i++ ) {
        itg += _c[i] * _term[i]->getItg();
    }
    return itg;
}

template <class T>
EvtValError EvtPdfSum<T>::compute_integral( int N ) const
{
    EvtValError itg( 0.0, 0.0 );
    for ( size_t i = 0; i < nTerms(); i++ )
        itg += _c[i] * _term[i]->getItg( N );
    return itg;
}

/*
 * Sample points randomly according to the sum of PDFs. First throw a random number uniformly
 * between zero and the value of the sum integral. Using this random number select one
 * of the PDFs. The generate a random point according to that PDF.
 */

template <class T>
T EvtPdfSum<T>::randomPoint()
{
    if ( !this->_itg.valueKnown() )
        this->_itg = compute_integral();

    double max = this->_itg.value();
    double rnd = EvtRandom::Flat( 0, max );

    double sum = 0.;
    size_t i;
    for ( i = 0; i < nTerms(); i++ ) {
        double itg = _term[i]->getItg().value();
        sum += _c[i] * itg;
        if ( sum > rnd )
            break;
    }

    return _term[i]->randomPoint();
}

#endif
