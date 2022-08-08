
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

#ifndef EVT_DALITZ_RES_PDF_HH
#define EVT_DALITZ_RES_PDF_HH

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtPdf.hh"

/*
 * Pole compensating function for terms that exibit a resonant structure
 * in one dimension only.
 *
 * f =  1    g*m0
 *     --  ------------------
 *     pi  (q-q0)^2 + g^2m0^2
 *
 * m is the mass of the resonance, g is its width. The approximation works well for a narrow
 * resonance. It is also readily integrable over the Dalitz plot coordinate to produce
 *
 * Int = 1/pi atan((q-q0)/(g*m0))
 */

class EvtDalitzResPdf : public EvtPdf<EvtDalitzPoint> {
  public:
    EvtDalitzResPdf( const EvtDalitzPlot& dp, double m0, double g0,
                     EvtCyclic3::Pair pairRes );

    EvtPdf<EvtDalitzPoint>* clone() const override
    {
        return new EvtDalitzResPdf( *this );
    }

    using EvtPdf<EvtDalitzPoint>::compute_integral;
    EvtValError compute_integral( int N ) const override;
    EvtDalitzPoint randomPoint() override;
    double pdfMaxValue() const;

  protected:
    double pdf( const EvtDalitzPoint& ) const override;

  private:
    EvtDalitzPlot _dp;
    double _m0;                // mass
    double _g0;                // width
    EvtCyclic3::Pair _pair;    // resonant pair
};

#endif
