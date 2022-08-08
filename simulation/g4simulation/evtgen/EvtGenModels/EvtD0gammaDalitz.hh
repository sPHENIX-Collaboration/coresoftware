
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

#ifndef __EVTD0GAMMADALITZ_HH__
#define __EVTD0GAMMADALITZ_HH__

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtSpinType.hh"

#include <vector>

class EvtParticle;

class EvtD0gammaDalitz : public EvtDecayAmp {
  private:
    int _d1;
    int _d2;
    int _d3;

    bool _isKsPiPi;

    // Useful constants.
    static const EvtSpinType::spintype& _SCALAR;
    static const EvtSpinType::spintype& _VECTOR;
    static const EvtSpinType::spintype& _TENSOR;

    static const EvtDalitzReso::CouplingType& _EtaPic;
    static const EvtDalitzReso::CouplingType& _PicPicKK;

    static const EvtDalitzReso::NumType& _RBW;
    static const EvtDalitzReso::NumType& _GS;
    static const EvtDalitzReso::NumType& _KMAT;

    static const EvtCyclic3::Pair& _AB;
    static const EvtCyclic3::Pair& _AC;
    static const EvtCyclic3::Pair& _BC;

    // Values to be read or computed based on values in the evt.pdl file.
    // IDs of the relevant particles.
    EvtId _BP;
    EvtId _BM;
    EvtId _B0;
    EvtId _B0B;
    EvtId _D0;
    EvtId _D0B;
    EvtId _KM;
    EvtId _KP;
    EvtId _K0;
    EvtId _K0B;
    EvtId _KL;
    EvtId _KS;
    EvtId _PIM;
    EvtId _PIP;

    // Flavor of the B mother.
    EvtId _bFlavor;

    // Masses of the relevant particles.
    double _mD0;
    double _mKs;
    double _mPi;
    double _mK;

    void readPDGValues();
    void reportInvalidAndExit() const;

    EvtComplex dalitzKsPiPi( const EvtDalitzPoint& point ) const;
    EvtComplex dalitzKsKK( const EvtDalitzPoint& point ) const;

  public:
    std::string getName() override;
    EvtDecayBase* clone() override;

    void init() override;
    void initProbMax() override;

    void decay( EvtParticle* p ) override;
};

#endif
