
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

#ifndef EVTVECTORISR_HH
#define EVTVECTORISR_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

class EvtParticle;

// Description:
//   This is a special decay model to generate e+e- -> phi gamma + soft gammas
//   using soft collinear ISR calculation from AfkQed
//   This is implemented as a decay of the VPHO.

class EvtVectorIsr : public EvtDecayIncoherent {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;

    void init() override;

    void initProbMax() override;

    double ckhrad1( double xx, double a, double b );

    void ckhrad( const double& e_beam, const double& q2_min, double& e01,
                 double& e02, double& f );

  private:
    double csfrmn, csbkmn;
    double fmax;
    bool firstorder;
};

#endif
