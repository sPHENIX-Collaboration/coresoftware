
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

#ifndef EVTLB2LLL_HH
#define EVTLB2LLL_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtTensor4C.hh"

#include "EvtGenModels/EvtWilsonCoefficients.hh"

// Desription: Routine to implement Lambda_b0 -> Lambda_0 l+ l- decays accroding to
//             several models: Chen. Geng.
//                             Aliev. Ozpineci. Savci.

class EvtLb2Lll : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;
    void init() override;
    void initProbMax() override;
    void calcAmp( EvtAmp* amp, EvtParticle* parent );

    EvtTensor4C EvtLeptonTG5Current( const EvtDiracSpinor& d,
                                     const EvtDiracSpinor& dp );

  private:
    double m_polarizationLambdab0;
    double m_maxProbability;
    double m_poleSize;
    long m_noTries;
    double m_omega;

    std::string m_decayName;
    std::string m_polarizationIntroduction;
    std::string m_HEPmodel;
    std::string m_FFtype;
    std::string m_effectContribution;

    EvtWilsonCoefficients m_WC;
};

#endif
