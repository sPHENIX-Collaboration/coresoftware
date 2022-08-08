
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

#ifndef EVT_PTO3P_AMP_FACTORY_HH
#define EVT_PTO3P_AMP_FACTORY_HH

#include "EvtGenBase/EvtAmpFactory.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"

#include <string>
#include <vector>

class EvtPto3PAmpFactory final : public EvtAmpFactory<EvtDalitzPoint> {
  public:
    EvtPto3PAmpFactory( const EvtDalitzPlot& dp ) :
        EvtAmpFactory<EvtDalitzPoint>(), _dp( dp )
    {
    }
    EvtPto3PAmpFactory( EvtPto3PAmpFactory&& ) = default;
    EvtPto3PAmpFactory( const EvtPto3PAmpFactory& ) = default;

    EvtAmpFactory<EvtDalitzPoint>* clone() const override
    {
        return new EvtPto3PAmpFactory( *this );
    }

    void processAmp( EvtComplex c, std::vector<std::string> vv,
                     bool conj ) override;

  private:
    double matchIsobarCoef( EvtAmplitude<EvtDalitzPoint>& amp,
                            EvtPdf<EvtDalitzPoint>& pdf, EvtCyclic3::Pair i );

    EvtDalitzPlot _dp;
};

#endif
