
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

#ifndef EVTDDALITZ_HH
#define EVTDDALITZ_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtFlatte.hh"

#include <vector>

class EvtParticle;

class EvtDDalitz : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;

  private:
    int _d1, _d2, _d3, _flag;

    EvtComplex amplDtoK0PiPi( EvtVector4R p4_p, EvtVector4R moms1,
                              EvtVector4R moms2, EvtVector4R moms3 );
    EvtComplex amplDtoK0KK( EvtVector4R p4_p, EvtVector4R moms1,
                            EvtVector4R moms2, EvtVector4R moms3 );

    vector<EvtFlatteParam> _kkpi_params;
};

#endif
