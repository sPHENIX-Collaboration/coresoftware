
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

#ifndef EVT_DALITZ_FLAT_PDF_HH
#define EVT_DALITZ_FLAT_PDF_HH

#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtPdf.hh"

#include <assert.h>

/*
 * Uniform PDF defined on a Dalitz plot.
 */

class EvtDalitzFlatPdf : public EvtPdf<EvtDalitzPoint> {
  public:
    EvtDalitzFlatPdf( const EvtDalitzPlot& dp );
    EvtDalitzFlatPdf( const EvtDalitzFlatPdf& other );
    EvtPdf<EvtDalitzPoint>* clone() const override;

    using EvtPdf<EvtDalitzPoint>::compute_integral;
    EvtValError compute_integral( int N ) const override;
    EvtDalitzPoint randomPoint() override;

  protected:
    double pdf( const EvtDalitzPoint& ) const override;

    EvtDalitzPlot _dp;
};

#endif
