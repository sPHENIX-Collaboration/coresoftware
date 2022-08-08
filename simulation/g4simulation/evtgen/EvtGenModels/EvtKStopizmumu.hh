
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

#ifndef EVTKSTOPIZMUMU_HH
#define EVTKSTOPIZMUMU_HH

#include "EvtGenBase/EvtDecayAmp.hh"

#include <string>

class EvtParticle;

// Description: Routine to implement KS -> pi0 mu mu; see JHEP08(1998)004.

class EvtKStopizmumu : public EvtDecayAmp {
  public:
    std::string getName() override { return "KS_PI0MUMU"; }

    EvtDecayBase* clone() override { return new EvtKStopizmumu; }

    void init() override;

    void initProbMax() override { setProbMax( 1.0e-10 ); }

    void decay( EvtParticle* p ) override;

    double F_z( const double& z, const double& rvsq );
    EvtComplex G_z( const double& z );
    double Wpol_z( const double& z, const double& as, const double& bs );
    EvtComplex chi_z( const double& z, const double& rpisq );
    EvtComplex Wpipi_z( const double& z, const double& alpha_s,
                        const double& beta_s, const double& rvsq,
                        const double& rpisq, const double& z0 );
};

#endif    //EVTKTOPIZMUMU_HH
