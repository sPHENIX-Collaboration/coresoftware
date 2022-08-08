
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

#ifndef EVTGENMODELS_EVTBTODDALITZCPK_HH
#define EVTGENMODELS_EVTBTODDALITZCPK_HH 1

// Include files
#include "EvtGenBase/EvtDecayAmp.hh"

/** @class EvtBToDDalitzCPK EvtBToDDalitzCPK.hh EvtGenModels/EvtBToDDalitzCPK.hh
 *  Decay Model for B->DK, (adds the possibility to use D0->Ks pi pi to
 *  find gamma with a Dalitz analysis
 *
 *  @author Patrick Robbe
 *  @date   2003-12-08
 */

class EvtBToDDalitzCPK : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtBToDDalitzCPK* clone() override;

    void decay( EvtParticle* p ) override;
    void init() override;

    void initProbMax() override;

  private:
    int _flag;
};
#endif    // EVTGENMODELS_EVTBTODDALITZCPK_HH
