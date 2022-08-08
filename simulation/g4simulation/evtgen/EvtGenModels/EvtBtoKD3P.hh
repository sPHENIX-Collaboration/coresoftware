
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

#ifndef EVT_BTOKD3P
#define EVT_BTOKD3P

class EvtParticle;
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

// Decay model that does the decay B+ -> K+ D , D -> 3 psudoscalars.
//
// The B- daughters specified in the decay file should be K-, D0, D0,
// where the first D0 is produced via b->c decay and the second via b->u.
// In reality, only one D daughter exists, so the first two
// daughters must be defined to decay to the same final state using
// the EvtPto3P model, but with CP-conjugate amplitudes.
//
// For a given point in the Pto3P Dalitz plot,
// the total amplitude is \propto [A1 + A2 r exp(i(phase))], where
//
// A1 & A2 are the amplitudes of the D0 and D0bar to decay into that
// Dalitz plot point,
//
// r is the (positive) ratio between the A(B->B0bar K) and A(B->D0 K)
// B decay amplitudes,
//
// phase is the total phase difference (weak phase + strong phase) between
// A(B->D0bar K) and A(B->B0 K).
//
// Note that this model knows nothing about your convention for the
// sign of the phase, so when specifying the decay of a B- you need to
// change the order of D0 and D0bar and change the total phase so that
// the sign of the weak phase flips with respect to the parameters of B+.

class EvtBtoKD3P : public EvtDecayAmp {
  public:
    EvtDecayBase* clone() override;

    // Initialize model
    void init() override;
    void initProbMax() override;
    void decay( EvtParticle* p ) override;

    // we really have two daughters, although three are listed in the .dec file:
    int nRealDaughters() override { return 2; }

    std::string getName() override;

  protected:
    // parameters:
    double _r;
    EvtComplex _exp;

    // other:
    const EvtDecayBase* _model1 = nullptr;
    const EvtDecayBase* _model2 = nullptr;
    bool _decayedOnce = false;
};

#endif
