
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

#ifndef EVT_PTO3P_AMP_SMPRSL_HH
#define EVT_PTO3P_AMP_SMPRSL_HH

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtPto3PAmp.hh"

class EvtComplex;

class EvtPto3PAmpSmpResolution : public EvtPto3PAmp {
  public:
    EvtPto3PAmpSmpResolution( EvtDalitzPlot dp, EvtCyclic3::Pair pairAng,
                              EvtCyclic3::Pair pairRes,
                              EvtSpinType::spintype spin,
                              const EvtPropagator& prop, NumType typeN );

    EvtAmplitude<EvtDalitzPoint>* clone() const override
    {
        return new EvtPto3PAmpSmpResolution( *this );
    }

    EvtComplex evalPropagator( double m ) const override;

    void setResolution( double bias, double sigma )
    {
        _bias = bias;
        _sigma = sigma;
    }

  private:
    double _bias;
    double _sigma;
};

#endif
