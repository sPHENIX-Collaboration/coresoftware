
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

#ifndef __EVTDALITZRESO_HH__
#define __EVTDALITZRESO_HH__

#include "EvtGenBase/EvtBlattWeisskopf.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzPoint.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"

#include <map>
#include <string>
#include <vector>

using std::map;
using std::vector;

class EvtBlattWeisskopf;

class EvtDalitzReso final {
  public:
    // Numerator type
    enum NumType
    {
        NBW = 0,
        RBW_ZEMACH = 1,
        RBW_KUEHN = 2,
        RBW_CLEO = 3,
        RBW_ZEMACH2 = 4,
        GS_CLEO = 5,
        K_MATRIX = 6,
        RBW_CLEO_ZEMACH = 7,
        GS_CLEO_ZEMACH = 8,
        LASS = 9,
        K_MATRIX_I = 10,
        K_MATRIX_II = 11,
        GAUSS_CLEO = 12,
        GAUSS_CLEO_ZEMACH = 13,
        FLATTE = 14,
        NON_RES = 15,
        NON_RES_LIN = 16,
        NON_RES_EXP = 17
    };

    // Coupling type
    //  ChgPion : pi+ pi-
    //  NeuPion : pi0 pi0
    //  Pion    : 0.5*[(pi+ pi-) + (pi0 pi0)]
    //  ChgKaon : K+ K-
    //  NeuKaon : K0 K0
    //  Kaon    : 0.5*[(K+ K-) + (K0 K0)]
    //  EtaPion : eta pi0
    enum CouplingType
    {
        Undefined = 0,
        PicPic = 1,
        PizPiz,
        PiPi,
        KcKc,
        KzKz,
        KK,
        EtaPic,
        EtaPiz,
        PicPicKK,
        WA76
    };

    EvtDalitzReso() : _typeN( NON_RES ){};

    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes,
                   NumType typeN, double alpha = 0.0 ) :
        _dp( dp ), _pairRes( pairRes ), _typeN( typeN ), _alpha( alpha ){};

    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairAng,
                   EvtCyclic3::Pair pairRes, EvtSpinType::spintype spin,
                   double m0, double g0, NumType typeN, double f_b = 0.0,
                   double f_d = 1.5 );

    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairAng,
                   EvtCyclic3::Pair pairRes, EvtSpinType::spintype spin,
                   double m0, double g0, NumType typeN, double m0_mix,
                   double g0_mix, double delta_mix, EvtComplex amp_mix );

    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairAng,
                   EvtCyclic3::Pair pairRes, EvtSpinType::spintype spin,
                   double m0, NumType typeN, double g1, double g2,
                   CouplingType coupling2 );

    // K-matrix
    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes,
                   std::string nameIndex, NumType typeN, EvtComplex fr12prod,
                   EvtComplex fr13prod, EvtComplex fr14prod,
                   EvtComplex fr15prod, double s0prod );

    // LASS
    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes, double m0,
                   double g0, double a, double r, double B, double phiB, double R,
                   double phiR, double cutoff = -1, bool scaleByMOverQ = false );

    //Flatte
    EvtDalitzReso( const EvtDalitzPlot& dp, EvtCyclic3::Pair pairRes, double m0 );

    EvtDalitzReso* clone() const { return new EvtDalitzReso( *this ); }

    EvtComplex evaluate( const EvtDalitzPoint& p );

    void set_fd( double R ) { _vd.set_f( R ); }
    void set_fb( double R ) { _vb.set_f( R ); }

    void addFlatteParam( const EvtFlatteParam& param )
    {
        _flatteParams.push_back( param );
    }

  private:
    EvtComplex psFactor( double& ma, double& mb, double& m );
    EvtComplex psFactor( double& ma1, double& mb1, double& ma2, double& mb2,
                         double& m );
    EvtComplex propGauss( const double& m0, const double& s0, const double& m );
    EvtComplex propBreitWigner( const double& m0, const double& g0,
                                const double& m );
    EvtComplex propBreitWignerRel( const double& m0, const double& g0,
                                   const double& m );
    EvtComplex propBreitWignerRel( const double& m0, const EvtComplex& g0,
                                   const double& m );
    EvtComplex propBreitWignerRelCoupled( const double& m0, const EvtComplex& g1,
                                          const EvtComplex& g2, const double& m );
    EvtComplex propGounarisSakurai( const double& m0, const double& g0,
                                    const double& k0, const double& m,
                                    const double& g, const double& k );
    inline double GS_f( const double& m0, const double& g0, const double& k0,
                        const double& m, const double& k );
    inline double GS_h( const double& m, const double& k );
    inline double GS_dhods( const double& m0, const double& k0 );
    inline double GS_d( const double& m0, const double& k0 );

    EvtComplex numerator( const EvtDalitzPoint& p, const EvtTwoBodyKine& vb,
                          const EvtTwoBodyKine& vd );
    double angDep( const EvtDalitzPoint& p );
    EvtComplex mixFactor( EvtComplex prop, EvtComplex prop_mix );
    EvtComplex Fvector( double s, int index );
    EvtComplex lass( double s );
    EvtComplex flatte( const double& m );

    inline EvtComplex sqrtCplx( double in )
    {
        return ( in > 0 ) ? EvtComplex( sqrt( in ), 0 )
                          : EvtComplex( 0, sqrt( -in ) );
    }

    // Dalitz plot
    EvtDalitzPlot _dp;

    // Pairing indices:
    EvtCyclic3::Pair _pairAng;    // angular
    EvtCyclic3::Pair _pairRes;    // resonance

    // Spin
    EvtSpinType::spintype _spin;

    // Numerator type
    NumType _typeN;

    // Nominal mass and width
    double _m0, _g0;

    // Vertices
    EvtTwoBodyVertex _vb;
    EvtTwoBodyVertex _vd;

    // Daughter masses
    double _massFirst, _massSecond;

    // variables for electromagnetic mass mixing
    double _m0_mix, _g0_mix, _delta_mix;
    EvtComplex _amp_mix;

    // variables for coupled Breit-Wigner
    double _g1, _g2;
    CouplingType _coupling2;

    // variables for Blatt-Weisskopf form factors
    double _f_b, _f_d;

    // K-matrix
    int _kmatrix_index;
    EvtComplex _fr12prod, _fr13prod, _fr14prod, _fr15prod;
    double _s0prod;

    // LASS
    double _a;
    double _r;
    double _Blass;
    double _phiB;
    double _R;
    double _phiR;
    double _cutoff;
    bool _scaleByMOverQ;

    //Nonresonant
    double _alpha;

    // Flatte
    std::vector<EvtFlatteParam> _flatteParams;
};

#endif
