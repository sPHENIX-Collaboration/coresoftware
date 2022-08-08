
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

#ifndef EvtBcVNpi_HH
#define EvtBcVNpi_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"

#include "EvtGenModels/EvtBCVFF.hh"
#include "EvtGenModels/EvtWnPi.hh"

#include <iostream>
#include <memory>

using std::endl;
using std::string;

class EvtBcVNpi : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;
    void initProbMax() override;
    void init() override;
    void decay( EvtParticle* p ) override;

  protected:
    int nCall;
    int whichfit, idVector;
    std::unique_ptr<EvtBCVFF> ffmodel;
    std::unique_ptr<EvtWnPi> wcurr;

    EvtComplex Fpi( EvtVector4R q1, EvtVector4R q2 );
    EvtComplex BWa( EvtVector4R q );
    EvtComplex BWf( EvtVector4R q );
    EvtComplex BWr( EvtVector4R q );
    EvtVector4C JB( EvtVector4R q1, EvtVector4R q2, EvtVector4R q3,
                    EvtVector4R q4, EvtVector4R q5 );
};
#endif
