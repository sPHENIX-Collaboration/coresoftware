
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

#ifndef EVTLAMBDAPBARGAMMA_HH
#define EVTLAMBDAPBARGAMMA_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"

// Description:Implementation of the decay B- -> lambda p_bar gamma according to
// Cheng, Yang; hep-ph/0201015

class EvtLambdaP_BarGamma : public EvtDecayAmp {
  public:
    EvtLambdaP_BarGamma();
    ~EvtLambdaP_BarGamma() { ; }

    std::string getName() override;
    EvtDecayBase* clone() override;
    void decay( EvtParticle* p ) override;
    void init() override;
    void initProbMax() override;

  private:
    // some constants to make the code easier to read and maintain
    // these three should be constants... (implementation of getMass() prohibits this)
    double _mLambdab;    //  = 5.624;                          // Lambda_b mass
    double _mLambda0;    //  = 1.115684;                       // Lambda0 mass
    double _c7Eff;    //  = -0.31;                          // Wilson coefficient
    double _mb;       //  =  4.4;                           // running b mass
    double _mV;    //  =  5.42;                          // pole mass vector current
    double _mA;    //  =  5.86;                          // pole mass axial current
    double _GF;    //  =  1.166E-5;                      // Fermi constant
    double _gLambdab;    // =  16;                            // coupling constant Lambda_b -> B- p
    double _e0;    //  =  1;                             // electromagnetic coupling (+1)
    double _g1;    //  =  0.64;                          // heavy-light form factors at q_mSqare
    double _g2;            //  = -0.10;
    double _f1;            //  =  0.64;
    double _f2;            //  = -0.31;
    double _VtbVtsStar;    // = 0.038;                          // |V_tb V_ts^*|

    // user never needs to call this -> private
    // baryonic form factors f(p), g(p), at p=0
    double f0( const double f_qm,
               int n = 1 ) const;    // calculate f(0) with f(q_max)
    double g0( const double f_qm,
               int n = 1 ) const;    // calculate g(0) with g(q_max)

    // shorthand for constants a and b in the formula
    double constA() const;
    double constB() const;

    // initialize phasespace and calculate the amplitude for one (i=0,1) state of the photon
    EvtComplex calcAmpliude( const EvtParticle* p, const unsigned int polState );
};

#endif
