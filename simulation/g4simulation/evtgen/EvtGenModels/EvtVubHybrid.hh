
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

#ifndef EVTVUBHYBRID_HH
#define EVTVUBHYBRID_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

#include "EvtGenModels/EvtVubdGamma.hh"

#include <memory>
#include <vector>

class EvtParticle;
class RandGeneral;

// Description:
// Class to generate inclusive B to X_u l nu decays.
// This class is based on EvtVub by Sven Menke with an update to
// generate the inclusive decays in such a way that the right
// mix of inclusive and exclusive decays is obtained:
// "Hybrid Model" by Dominique Fortin.
// NOTE:
// - A set of weights (for bins in the kinematic variables mX, q2, El)
//   is read from DECAY.DEC. This set of weights must be consistent
//   with the other parameters specified (excl. BF, non-res BF, mb, a).
// - If no binning/weights are specified in DECAY.DEC the hybrid
//   reweighting is not activated

class EvtVubHybrid : public EvtDecayIncoherent {
  public:
    std::string getName() override;

    EvtDecayBase* clone() override;

    void initProbMax() override;

    void init() override;

    void decay( EvtParticle* p ) override;

    void readWeights( int startArg = 0 );

    double getWeight( double mX, double q2, double El );

  private:
    double findPFermi();

    enum
    {
        nParameters = 3,
        nVariables = 3
    };

    bool _noHybrid =
        false;    // _noHybrid will be set TRUE if the DECAY.DEC file has no binning or weights
    bool _storeQplus =
        true;    // _storeQplus should alwasy be TRUE: writes out Fermi motion parameter

    double _mb = 4.62;        // the b-quark pole mass in GeV (try 4.65 to 4.9)
    double _a = 2.27;         // Parameter for the Fermi Motion (1.29 is good)
    double _alphas = 0.22;    // Strong Coupling at m_b (around 0.24)
    double _dGMax = 3.;       // max dGamma*p2 value;
    int _nbins = 0;
    double _masscut = 0.28;
    std::vector<double> _bins_mX;
    std::vector<double> _bins_q2;
    std::vector<double> _bins_El;
    std::vector<double> _weights;
    std::unique_ptr<EvtVubdGamma> _dGamma;    // calculates the decay rate
    std::vector<double> _pf;
};

#endif
