
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

#ifndef EVT_LASS_AMP_HH
#define EVT_LASS_AMP_HH

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"

#include <string>

class EvtComplex;

class EvtLASSAmp : public EvtAmplitude<EvtDalitzPoint> {
  public:
    EvtLASSAmp( EvtDalitzPlot* dp, EvtCyclic3::Pair pair, double m0, double g0,
                double a, double r, double cutoff, std::string subtype = "LASS" );

    EvtComplex amplitude( const EvtDalitzPoint& p ) const override;

    EvtAmplitude<EvtDalitzPoint>* clone() const override
    {
        return new EvtLASSAmp( *this );
    }

  private:
    EvtDalitzPlot* _dalitzSpace;

    EvtCyclic3::Pair _pair;

    double _m0;
    double _g0;
    double _q0;
    double _r;
    double _a;
    double _cutoff;
    std::string _subtype;
};

#endif
