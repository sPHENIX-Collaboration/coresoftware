
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

#ifndef EVTPI0DALITZ_HH
#define EVTPI0DALITZ_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtPi0Dalitz : public EvtDecayProb {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;

  private:
    double m_poleSize{ 0.00000002 };

    // Following are rho mass and width, but in order to keep consistency
    // with what was done before do not use data from particle table.
    const double m_m0Sq{ 0.768 * 0.768 };
    const double m_m0SqG0Sq{ m_m0Sq * 0.151 * 0.151 };
};

#endif
