
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

#ifndef EVTBTODIBARYONLNUPQCDFF_HH
#define EVTBTODIBARYONLNUPQCDFF_HH

class EvtParticle;
#include <vector>

// Description: Module for computation of B->ppbar form factors according
// to pQCD counting rules, see arXiv:1107.0801.

class EvtBToDiBaryonlnupQCDFF {
  public:
    class FormFactors final {
      public:
        double F1, F2, F3, F4, F5;
        double G1, G2, G3, G4, G5;
    };

    EvtBToDiBaryonlnupQCDFF();

    EvtBToDiBaryonlnupQCDFF( std::vector<double>& DParameters );

    void getDiracFF( EvtParticle* parent, double dibaryonMass,
                     EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const;

    void getRaritaFF( EvtParticle* parent, double dibaryonMass,
                      EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const;

    void getFF( EvtParticle* parent, double dibaryonMass,
                EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const;

  private:
    std::vector<double> DPars;
    int nDPars;
};

#endif
