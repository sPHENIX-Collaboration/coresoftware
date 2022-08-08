
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

#ifndef EVT_PTO3P_HH
#define EVT_PTO3P_HH

#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtIntervalDecayAmp.hh"

#include <vector>

class EvtPto3P : public EvtIntervalDecayAmp<EvtDalitzPoint> {
  public:
    EvtPto3P() {}
    ~EvtPto3P() {}
    std::string getName() override { return "PTO3P"; }
    EvtDecayBase* clone() override { return new EvtPto3P(); }

    EvtAmpFactory<EvtDalitzPoint>* createFactory(
        const EvtMultiChannelParser& parser ) override;
    std::vector<EvtVector4R> initDaughters( const EvtDalitzPoint& p ) const override;

    EvtDalitzPlot dp();
};

#endif
