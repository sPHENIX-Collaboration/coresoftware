
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

#ifndef EVTSVVNONCPEIGEN_HH
#define EVTSVVNONCPEIGEN_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

// Description: Routine to decay scalar -> vector vector
//              and has CP violation.
//
//              This model does all the ckm-suppressed decays and mixing for you. It randomly 'overwrites'
//              any reco or tagging state as set in the Y(4S) decay model (VSS_(B)MIX) with its own generated states.
//
//              As such, the corresponding dec file requires only one decay-mode description, for example:
//              Decay MyB0
//              1.000    rho+ MyD*-       SVV_NONCPEIGEN dm beta gamma 0.322 0.31 0.941 0 0.107 1.42 0.02 0 0.02 0 0.02 0 ;
//              EndDecay
//              and furthermore Y(4S) only needs to decay to B0's (or B0bar's).
//              The decay above should be a CKM-favored mode (eg. B0->D*-rho+ or B0bar->D*+rho-).
//              All ckm-suppressed decays and the mixing are derived from this line in the ::Decay function.
//
//              There are 15 or 27 arguments. The first three are dm, phase1
//              and phase2. dm is the B0-B0bar mass difference. Phases 1
//              and 2 are the CKM weak phases relevant for the particular mode,
//              eg for B-->DstRho phase1 is beta and phase2 is gamma.
//
//              The next arguments are the 2 amplitudes (= 12 input parameters)
//              in the order: A_f, Abar_f. In the example above, the 'A_f' amplitude now
//              stands for the ckm-favored decay 'B0->D*-rho+', and 'Abar_f' stands for 'B0bar->D*-rho+'
//
//              Each amplitude has its 3 helicity states in the order +, 0, -, which are each
//              specified by a magnitude and a strong phase.
//
//              The last 2 arguments A_fbar and Abar_fbar (=12 input parameters) are not necessary,
//              but can included if one wants to set them differently from A_f, Abar_f.
//
//              Mind you that Hbar_+- = H_-+ (ignoring the weak phase, which flips sign).
//              It is custumary to select one set of helicity states (eg H_+-) and to adopt these for
//              the CP-conjugate decays as well (ie. depict Hbar_-+ with H_+-), which is the interpretation
//              we use for the input-parameters above.
//              However, the angular decay in EvtGen is just a formula in which helicity amplitudes are 'plugged' in,
//              making no difference between B0 or B0bar decays. In the model below we (thus) account for the +-
//              flipping between B0 and B0bar.

class EvtSVVNONCPEIGEN : public EvtDecayAmp {
  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void initProbMax() override;
    void init() override;

    void decay( EvtParticle* p ) override;

    std::string getParamName( int i ) override;
    std::string getParamDefault( int i ) override;

  private:
    EvtComplex _A_f[12];
};

#endif
