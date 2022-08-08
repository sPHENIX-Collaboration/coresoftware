
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

#ifndef EVTPSI2JPSIPIPI_HH
#define EVTPSI2JPSIPIPI_HH

#include "EvtGenBase/EvtDecayAmp.hh"

#include <array>
#include <string>

class EvtDecayBase;
class EvtParticle;

// Description: Header file for the model "PSI2JPSIPIPI" which generates
//              psi2S -> J/psi pi+ pi- decays based on hep-ph/1507.07985

class EvtPsi2JpsiPiPi : public EvtDecayAmp {
  public:
    EvtPsi2JpsiPiPi();

    std::string getName() override;
    EvtDecayBase* clone() override;
    void initProbMax() override;
    void init() override;
    void decay( EvtParticle* p ) override;

  private:
    bool tree;
    double phi;    // LO vs NLO mixing angle (radians)
    double cosPhi, cos2Phi, sinPhi, sin2Phi;
    // NLO corrections
    static const int nQ = 6;    // number of terms in mPiPi interpolation
    std::array<double, nQ> c0, c1, c2, s1, s2;

    void setNLOArrays();
};

#endif
