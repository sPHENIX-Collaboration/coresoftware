
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

#ifndef EVTSLBARYONAMP_HH
#define EVTSLBARYONAMP_HH

#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtAmp;
class EvtComplex;
class EvtDiracSpinor;
class EvtParticle;
class EvtRaritaSchwinger;
class EvtSemiLeptonicFF;
class EvtVector4C;
class EvtVector4R;

class EvtSLBaryonAmp : public EvtSemiLeptonicAmp {
  public:
    ~EvtSLBaryonAmp();

    //Daughters are initialized and have been added to the parent.
    //No need to carry around the daughters seperately!
    void CalcAmp( EvtParticle* parent, EvtAmp& amp,
                  EvtSemiLeptonicFF* FormFactors ) override;

    void CalcAmp( EvtParticle* parent, EvtAmp& amp,
                  EvtSemiLeptonicFF* FormFactors, EvtComplex r00,
                  EvtComplex r01, EvtComplex r10, EvtComplex r11 );

    double CalcMaxProb( EvtId parent, EvtId meson, EvtId lepton, EvtId nudaug,
                        EvtSemiLeptonicFF* FormFactors, EvtComplex r00,
                        EvtComplex r01, EvtComplex r10, EvtComplex r11 );

  private:
    EvtVector4C EvtBaryonVACurrent( const EvtDiracSpinor& Bf,
                                    const EvtDiracSpinor& Bi,
                                    EvtVector4R parent, EvtVector4R daught,
                                    const double* ff, int pflag );

    EvtVector4C EvtBaryonVARaritaCurrent( const EvtRaritaSchwinger& Bf_vect,
                                          const EvtDiracSpinor& Bi,
                                          EvtVector4R parent, EvtVector4R daught,
                                          const double* ff, int pflag );
};

#endif
