
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

#ifndef EVT_NONRESONANT_AMP_HH
#define EVT_NONRESONANT_AMP_HH

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtPto3PAmp.hh"
#include "EvtGenBase/EvtSpinType.hh"

class EvtComplex;

class EvtNonresonantAmp : public EvtAmplitude<EvtDalitzPoint> {
  public:
    EvtNonresonantAmp( EvtDalitzPlot* dp, EvtPto3PAmp::NumType type,
                       EvtCyclic3::Pair pair1, double par1 = 0,
                       EvtCyclic3::Pair pair2 = EvtCyclic3::AB, double par2 = 0,
                       EvtSpinType::spintype spin = EvtSpinType::SCALAR );

    EvtComplex amplitude( const EvtDalitzPoint& p ) const override;

    EvtAmplitude<EvtDalitzPoint>* clone() const override
    {
        return new EvtNonresonantAmp( *this );
    }

  private:
    EvtDalitzPlot* _dalitzSpace;

    EvtPto3PAmp::NumType _type;

    EvtCyclic3::Pair _pair1, _pair2;

    double _par1, _par2;

    EvtSpinType::spintype _spin;
};

#endif
