
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

#ifndef __EVTD0MIXDALITZ_HH__
#define __EVTD0MIXDALITZ_HH__

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinType.hh"

// Description:
//   The D0mixDalitz model, with many resonances and mixing implemented.

class EvtD0mixDalitz : public EvtDecayAmp {
  private:
    int _d1;
    int _d2;
    int _d3;

    // Mixing parameters.
    double _x;
    double _y;

    // q/p CP violation in the mixing.
    EvtComplex _qp;

    // Checker of the decay mode.
    bool _isKsPiPi;
    bool _isRBWmodel;

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

    // Masses of the relevant particles.
    double _mD0;
    double _mKs;
    double _mPi;
    double _mK;

    // Life time and decay rate.
    double _ctau;
    double _gamma;

    // Some useful integrals over the Dalitz plot.
    EvtComplex _iChi;
    EvtComplex _iChi2;

    void readPDGValues();
    EvtComplex dalitzKsPiPi( const EvtDalitzPoint& point );
    EvtComplex dalitzKsKK( const EvtDalitzPoint& point );

    // Time evolution functions for hamiltonian eigenstates.
    //    Negative exponential part removed.
    EvtComplex h1( const double& ct ) const;
    EvtComplex h2( const double& ct ) const;

    void reportInvalidAndExit() const
    {
        EvtGenReport( EVTGEN_ERROR, "EvtD0mixDalitz" )
            << "EvtD0mixDalitz: Invalid mode." << std::endl;
        exit( 1 );
    }

  public:
    EvtD0mixDalitz() :
        _d1( 0 ),
        _d2( 0 ),
        _d3( 0 ),
        _x( 0. ),
        _y( 0. ),
        _qp( 1. ),
        _isKsPiPi( false ),
        _isRBWmodel( true )
    {
    }

    // One-line inline functions.
    std::string getName() override { return "D0MIXDALITZ"; }
    EvtDecayBase* clone() override { return new EvtD0mixDalitz; }
    void initProbMax() override { setProbMax( 5200. ); }

    void init() override;
    void decay( EvtParticle* p ) override;
};

#endif
