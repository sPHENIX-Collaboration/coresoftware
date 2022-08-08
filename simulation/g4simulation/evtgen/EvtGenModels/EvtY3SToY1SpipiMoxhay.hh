
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

#ifndef EVTY3STOY1SPIPIMOXHAY_HH
#define EVTY3STOY1SPIPIMOXHAY_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

// Description: This model is based on the proposal by Tuan and Lipkin
//              (Phys.Lett.B206:349-353,1988) and the subsequent model
//              by Moxhay (Phys.Rev.D39:3497,1989) for the dipion spectrum
//              in Y(3S) -> pi+ pi- Y(1S). Please Note: in Moxhay's paper,
//              he wrote the fitted value of the parameter Im(B)/A as
//              -0.2983. However, using his quoted value leads to the wrong
//              spectrum. Changing the sign of his quoted Im(B)/A fixes the
//              shape and reproduces his result. Therefore, please pass
//              Im(B)/A = 0.2983 and Re(B)/A = 0.2196 to get the correct shape
//              based on his fit to the CLEO data.
//
// Example:
//
// Decay  Upsilon(3S)
//  1.0000    Upsilon  pi+  pi-     Y3STOY1SPIPIMOXHAY 0.2196 0.2983;
// Enddecay
//
//   --> the order of parameters is: Re(B)/A Im(B)/A

class EvtY3SToY1SpipiMoxhay : public EvtDecayProb {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void init() override;
    void initProbMax() override;
};

#endif
