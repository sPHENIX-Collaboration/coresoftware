
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

#ifndef EVTTAUHADNUKS_HH
#define EVTTAUHADNUKS_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtTauHadnu : public EvtDecayAmp {
  public:
    EvtTauHadnu() {}

    std::string getName() override;
    EvtDecayBase* clone() override;

    void initProbMax() override;
    void init() override;
    void decay( EvtParticle* p ) override;

  private:
    double _beta;
    double _mRho;
    double _gammaRho;
    double _mRhopr;
    double _gammaRhopr;
    double _mA1;
    double _gammaA1;

    double gFunc( double m2, int dupD );
    EvtComplex Fpi( double s, double xm1, double xm2 );
    EvtComplex BW( double s, double m, double gamma, double xm1, double xm2 );
};

#endif
