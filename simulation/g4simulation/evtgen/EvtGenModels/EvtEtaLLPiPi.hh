
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

#ifndef EVT_ETALLPIPI_HH
#define EVT_ETALLPIPI_HH

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDecayProb.hh"

#include <string>

class EvtParticle;

// eta' -> mu+ mu- pi+ pi- or e+ e- pi+ pi-
// From Zhang Zhen-Yu et al, Chinese Phys. C 36, p926, 2012

class EvtEtaLLPiPi : public EvtDecayProb {
  public:
    EvtEtaLLPiPi() = default;

    void init() override;
    void initProbMax() override;

    std::string getName() override;
    EvtDecayBase* clone() override;

    void decay( EvtParticle* p ) override;

  private:
    void updateMassPars( double mLep, double mPi );

    double rhoWidth( double s, double m ) const;

    double F0( double sLL, double sPiPi ) const;

    double lambda( double a, double b, double c ) const;

    double ampSquared( EvtParticle* p ) const;

    double m_alpha{ 1.0 / 137.0 };
    double m_eSq{ 4.0 * EvtConst::pi * m_alpha };
    double m_fPi{ 0.0924 };
    double m_f8{ 1.3 * m_fPi };
    double m_f0{ 1.04 * m_fPi };
    double m_thetaMix{ 20.0 * EvtConst::pi / 180.0 };
    double m_mixSq{ 0.0 };
    double m_c1{ 1.0 };
    double m_c2{ 0.0 };
    double m_c3{ m_c1 - m_c2 };    // Eq 9
    double m_par1{ 1.0 - ( 3.0 * ( m_c1 - m_c2 + m_c3 ) / 4.0 ) };
    double m_parLL{ 3.0 * ( m_c1 - m_c2 - m_c3 ) / 4.0 };
    double m_parPiPi{ 3.0 * m_c3 / 2.0 };
    double m_rhoMass{ 0.775 };    // updated in init()
    double m_rhoMassSq{ m_rhoMass * m_rhoMass };
    double m_rhoGamma{ 0.149 };    // updated in init()
    double m_lepMass{ 0.106 };     // modified in updateMassPars()
    double m_lepMassSq{ m_lepMass * m_lepMass };
    double m_piMass{ 0.140 };    // modified in updateMassPars()
    double m_piMassSq{ m_piMass * m_piMass };
    double m_4LepMassSq{ 4.0 * m_lepMassSq };
    double m_4PiMassSq{ 4.0 * m_piMassSq };
};

#endif
