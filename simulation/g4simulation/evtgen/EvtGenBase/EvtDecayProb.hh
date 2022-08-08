
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

#ifndef EvtDecayProb_HH
#define EvtDecayProb_HH

#include "EvtGenBase/EvtDecayBase.hh"

class EvtParticle;

class EvtDecayProb : public EvtDecayBase {
  public:
    void makeDecay( EvtParticle* p, bool recursive = true ) override;

    void setProb( double prob ) { _prob = prob; }
    double getProb() { return _prob; }
    inline void setWeight( double weight ) { _weight = weight; }

    virtual ~EvtDecayProb() {}

  private:
    double _prob;
    double _weight;
};

#endif
