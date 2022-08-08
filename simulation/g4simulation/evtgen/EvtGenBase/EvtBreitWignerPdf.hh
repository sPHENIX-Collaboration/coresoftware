
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

#ifndef EVT_BREIT_WIGNER_PDF_HH
#define EVT_BREIT_WIGNER_PDF_HH

#include "EvtGenBase/EvtIntegPdf1D.hh"

// Breit-Wigner PDF

class EvtBreitWignerPdf : public EvtIntegPdf1D {
  public:
    EvtBreitWignerPdf( double min, double max, double m0, double g0 );
    EvtBreitWignerPdf( const EvtBreitWignerPdf& other );

    double pdf( const EvtPoint1D& x ) const override;
    EvtPdf<EvtPoint1D>* clone() const override
    {
        return new EvtBreitWignerPdf( *this );
    }

    double pdfIntegral( double m ) const override;
    double pdfIntegralInverse( double x ) const override;

    // accessors

    inline double m0() const { return _m0; }
    inline double g0() const { return _g0; }

  private:
    double _m0;
    double _g0;
};

#endif
