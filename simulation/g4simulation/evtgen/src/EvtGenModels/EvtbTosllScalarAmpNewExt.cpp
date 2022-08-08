
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

#include "EvtGenModels/EvtbTosllScalarAmpNewExt.hh"

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtbTosllAmpNewExt.hh"
#include "EvtGenModels/EvtbTosllFFNew.hh"
#include "EvtGenModels/EvtbTosllWilsCoeffNLO.hh"

#include <cstdlib>

//
// The main functiom for the amplitude calculation
//
// input:   *parent      - the pointer to the parent particle (B-meson, the
//                                          object of the EvtParticle class);
//          *formFactors - the pointer to instance of EvtbTosllFFNew class object;
//          *WilsCoeff   - the pointer to
//           mu          - the scale parameter, GeV;
//           Nf          - number of "effective" flavors (for b-quark Nf=5);
//           res_swch    - resonant switching parameter:
//                         = 0 the resonant contribution switched OFF,
//                         = 1 the resonant contribution switched ON;
//           ias         - switching parameter for \alpha_s(M_Z) value:
//                         = 0 PDG 1sigma minimal alpha_s(M_Z),
//                         = 1 PDG average value  alpha_s(M_Z),
//                         = 2 PDG 1sigma maximal alpha_s(M_Z).
//           Wolfenstein parameterization for CKM matrix
//                         CKM_A, CKM_lambda, CKM_barrho, CKM_bareta
//
// return:  amp          - amplitude for the decay  B -> P ell^+ ell^-
//
// Note: in our calculations we assume, that pseudoscalar meson is the first
//       daughter particle (iP=0) and leptons are the second and thirds
//       daughter particles (il1=1 and il2=2).
//
void EvtbTosllScalarAmpNewExt::CalcAmp(
    EvtParticle* parent, EvtAmp& amp, EvtbTosllFFNew* formFactors,
    EvtbTosllWilsCoeffNLO* WilsCoeff, double mu, int Nf, int res_swch, int ias,
    double CKM_A, double CKM_lambda, double CKM_barrho, double CKM_bareta,
    double ReA7, double ImA7, double ReA10, double ImA10 )
{
    //  FILE *mytest;

    EvtComplex unit1( 1.0, 0.0 );    // real unit
    EvtComplex uniti( 0.0, 1.0 );    // imaginary unit

    EvtComplex A7 = ReA7 * unit1 + ImA7 * uniti;
    EvtComplex A10 = ReA10 * unit1 + ImA10 * uniti;

    int iP = 0;    // pseudoscalar meson is the first daughter particle
    int il1 = 1,
        il2 = 2;    // leptons are the second and thirds daughter particles

    // transition momentum of the leptonic pair q=k1+k2 or q=p1-p2
    EvtVector4R q = parent->getDaug( il1 )->getP4() +
                    parent->getDaug( il2 )->getP4();

    // Mandelstam variable t=q^2
    double q2 = q.mass2();

    double M1 = parent->mass();                    // B - meson mass
    double M2 = parent->getDaug( iP )->mass();     // pseudoscalar meson mass
    double ml = parent->getDaug( il1 )->mass();    // leptonic mass
    double ms = 0.0;    // light quark mass from the dispersion QM
    double mc = formFactors->getQuarkMass( 4 );    // m_c mass from the dispersion QM
    double mb = formFactors->getQuarkMass( 5 );    // m_b mass from the dispersion QM
    //  double Mw = EvtPDL::getNominalMass("W+");  // W-boson mass
    //  double mt = EvtPDL::getNominalMass("t");   // t-quark mass
    double Mw = 80.403;    // GeV W-boson mass
    double mt = 174.2;     // GeV t-quark mass

    EvtComplex Vtb, Vtq, Vub, Vuq;    // V_{tb}, V_{tq}, V_{ub} and V_{uq}
    EvtComplex CKM_factor;            // V^*_{tq}*V_{tb}, where q={d,s}
    EvtComplex lambda_qu;    // V^*_{uq}*V_{ub}/V^*_{tq}*V_{tb}, where q={d,s}
    double Relambda_qu, Imlambda_qu;

    EvtId idparent = parent->getId();                   // B-meson Id
    EvtId iddaught = parent->getDaug( iP )->getId();    // The pseudoscalar meson Id

    // set of the light quark mass value
    if ( ( idparent == EvtPDL::getId( std::string( "B+" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K+" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B-" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K-" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "anti-K0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta'" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta'" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "f_0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "f_0" ) ) ) ) {
        ms = formFactors->getQuarkMass( 3 );    // m_s mass from the dispersion QM
        // V_{ts}
        Vtq = unit1 * ( 1.0 - 0.5 * pow( CKM_lambda, 2.0 ) ) +
              pow( CKM_lambda, 2.0 ) *
                  ( CKM_barrho * unit1 + CKM_bareta * uniti ) /
                  sqrt( 1.0 - pow( CKM_lambda, 2.0 ) );
        Vtq = -CKM_A * pow( CKM_lambda, 2.0 ) * Vtq;
        // V_{us}
        Vuq = CKM_lambda * unit1;
    }

    if ( ( idparent == EvtPDL::getId( std::string( "B+" ) ) &&
           iddaught == EvtPDL::getId( std::string( "pi+" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B-" ) ) &&
           iddaught == EvtPDL::getId( std::string( "pi-" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "pi0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "pi0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta'" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "eta'" ) ) ) ) {
        ms = formFactors->getQuarkMass( 2 );    // m_d mass from the dispersion QM
        // V_{td}
        Vtq = unit1 - ( 1.0 - 0.5 * pow( CKM_lambda, 2.0 ) ) *
                          ( CKM_barrho * unit1 + CKM_bareta * uniti ) /
                          sqrt( 1.0 - pow( CKM_lambda, 2.0 ) );
        Vtq = CKM_A * pow( CKM_lambda, 3.0 ) * Vtq;
        // V_{ud}
        Vuq = unit1 * ( 1.0 - 0.5 * pow( CKM_lambda, 2.0 ) -
                        0.125 * pow( CKM_lambda, 4.0 ) );
    }

    if ( ms < 0.001 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n The function EvtbTosllScalarAmpNew::CalcAmp(...)"
            << "\n Error in the model set!"
            << " ms = " << ms << std::endl;
        ::abort();
    }

    Vtb = unit1 * ( 1.0 - 0.5 * pow( CKM_A * CKM_lambda * CKM_lambda,
                                     2.0 ) );    // V_{tb}
    Vub = CKM_A * pow( CKM_lambda, 3.0 ) *
          ( CKM_barrho * unit1 - CKM_bareta * uniti ) /
          sqrt( 1.0 - pow( CKM_lambda, 2.0 ) );    // V_{ub}

    CKM_factor = conj( Vtq ) * Vtb;    // V^*_{tq}*V_{tb}

    lambda_qu = conj( Vuq ) * Vub /
                CKM_factor;    // V^*_{uq}*V_{ub}/V^*_{tq}*V_{tb}
    Relambda_qu = real( lambda_qu );
    Imlambda_qu = imag( lambda_qu );

    double fp, f0, ft;    // B -> P transition form-factors

    // To get the B -> P transition form-factors
    formFactors->getScalarFF( parent->getId(), parent->getDaug( iP )->getId(),
                              q2, fp, f0, ft );

    // The Wilson Coefficients preparation according to the paper
    // A.J.Buras, M.Munz, Phys.Rev.D52, p.189 (1995)
    EvtComplex c7gam = WilsCoeff->GetC7Eff( mu, Mw, mt, Nf, ias );
    c7gam = c7gam * A7;
    EvtComplex c9eff_b2q = WilsCoeff->GetC9Eff( 0, res_swch, ias, Nf, q2, mb,
                                                ms, mc, mu, mt, Mw, ml,
                                                Relambda_qu, Imlambda_qu );
    EvtComplex c9eff_barb2barq = WilsCoeff->GetC9Eff( 1, res_swch, ias, Nf, q2,
                                                      mb, ms, mc, mu, mt, Mw, ml,
                                                      Relambda_qu, Imlambda_qu );
    EvtComplex c10a = WilsCoeff->GetC10Eff( mt, Mw );
    c10a = c10a * A10;

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << "\n\n The function EvtbTosllScalarAmpNew::CalcAmp(...) passed."
    //      << "\n Particle masses:"
    //      << "\n B - meson mass M1 = " << M1
    //      << "\n P - meson mass M2 = " << M2
    //      << "\n leptonic mass  ml = " << ml
    //      << "\n light quark mass  = " << ms
    //      << "\n c - quark mass mc = " << mc
    //      << "\n b - quark mass mb = " << mb
    //      << "\n t - quark mass mt = " << mt
    //      << "\n W - boson mass Mw = " << Mw
    //      << "\n ============================================================================"
    //      << "\n Input parameters:"
    //      << "\n scale parameter        mu = " << mu
    //      << "\n number of flavors      Nf = " << Nf
    //      << "\n resonant switching        = " << res_swch
    //      << "\n parameter for alpha_s(M_Z) = " << ias
    //      << "\n ============================================================================"
    //      << "\n Vector form-factors at q^2 = " << q2
    //      << "   for B -> P transition:"
    //      << "\n fp = " << fp
    //      << "\n f0 = " << f0
    //      << "\n ft = " << ft
    //      << "\n ============================================================================"
    //      << "\n Wilson Coefficients:"
    //      << "\n Re(c7gam) = " << real(c7gam) << " Im(c7gam) = " << imag(c7gam)
    //      << "\n Re(c9eff_b2q) = " << real(c9eff_b2q)
    //        << " Im(c9eff_b2q) = " << imag(c9eff_b2q)
    //      << "\n Re(c9eff_barb2barq) = " << real(c9eff_barb2barq)
    //        << " Im(c9eff_barb2barq) = " << imag(c9eff_barb2barq)
    //      << "\n Re(c10a)  = " << real(c10a)  << " Im(c10a)  = " << imag(c10a)
    //      << std::endl;

    //      mytest = fopen("scalaroutput.txt","a");
    //      if(mytest != NULL){
    //	fprintf(mytest,"%lf\n",q2);
    //	fclose(mytest);
    //      }
    //      else{
    //         EvtGenReport(EVTGEN_ERROR,"EvtGen") << "\n Error in writing to file.\n"
    //         << std::endl;
    //	 return;
    //      }

    // 4- momentum of the B-meson in the the B-meson rest frame
    EvtVector4R p1 = parent->getP4Restframe();
    EvtVector4R hatp1 = p1 / M1;
    // 4-momentum of the pseudoscalar meson in the B-meson rest frame
    EvtVector4R p2 = parent->getDaug( 0 )->getP4();
    EvtVector4R hatp2 = p2 / M1;

    // 4-vector \hat q = q/M1
    EvtVector4R hatq = q / M1;
    // 4-vector \hat P= (p1 + p2)/M1
    EvtVector4R hatP = hatp1 + hatp2;

    double hats = q2 / pow( M1, 2 );
    double hatM2 = M2 / M1;
    double hatmb = mb / M1;
    double hatms = ms / M1;

    // Hadronic matrix element with m_s.NE.0
    EvtComplex a_b2q, a_barb2barq, b_b2q, b_barb2barq, c, d;

    a_b2q = c9eff_b2q * fp -
            2.0 * c7gam * ( hatmb + hatms ) * ft / ( 1.0 + hatM2 );
    a_barb2barq = c9eff_barb2barq * fp -
                  2.0 * c7gam * ( hatmb + hatms ) * ft / ( 1.0 + hatM2 );

    b_b2q = ( c9eff_b2q * ( f0 - fp ) +
              2.0 * c7gam * ( hatmb + hatms ) * ft / ( 1.0 + hatM2 ) ) *
            ( 1 - pow( hatM2, 2.0 ) ) / hats;
    b_barb2barq = ( c9eff_barb2barq * ( f0 - fp ) +
                    2.0 * c7gam * ( hatmb + hatms ) * ft / ( 1.0 + hatM2 ) ) *
                  ( 1 - pow( hatM2, 2.0 ) ) / hats;

    c = c10a * fp;

    d = c10a * ( 1.0 - pow( hatM2, 2 ) ) * ( f0 - fp ) / hats;

    // to find ell^+ and ell^- in the B-meson  daughters
    int charge1 = EvtPDL::chg3( parent->getDaug( 1 )->getId() );
    int charge2 = EvtPDL::chg3( parent->getDaug( 2 )->getId() );

    EvtParticle* lepPlus = 0;
    EvtParticle* lepMinus = 0;

    lepPlus = ( charge1 > charge2 ) ? parent->getDaug( 1 ) : parent->getDaug( 2 );
    lepMinus = ( charge1 < charge2 ) ? parent->getDaug( 1 )
                                     : parent->getDaug( 2 );

    EvtVector4C T1, T2;    // hadronic matrix element vector structures

    EvtVector4C lvc11, lvc12;    // spin structures for
    EvtVector4C lvc21, lvc22;    // the leptonic vector current

    EvtVector4C lac11, lac12;    // spin structures for
    EvtVector4C lac21, lac22;    // the leptonic axial current

    // B - and barB - mesons descriptors
    EvtIdSet bmesons( "B-", "anti-B0", "anti-B_s0", "B_c-" );
    EvtIdSet bbarmesons( "B+", "B0", "B_s0", "B_c+" );

    EvtId parentID = parent->getId();

    if ( bmesons.contains( parentID ) ) {
        // The amplitude for the decay barB -> barP ell^+ ell^-
        // (b -> q ell^+ ell^- transition)

        T1 = a_b2q * hatP + b_b2q * hatq;

        T2 = c * hatP + d * hatq;

        lvc11 = EvtLeptonVCurrent( lepPlus->spParent( 0 ),
                                   lepMinus->spParent( 0 ) );
        lvc21 = EvtLeptonVCurrent( lepPlus->spParent( 1 ),
                                   lepMinus->spParent( 0 ) );
        lvc12 = EvtLeptonVCurrent( lepPlus->spParent( 0 ),
                                   lepMinus->spParent( 1 ) );
        lvc22 = EvtLeptonVCurrent( lepPlus->spParent( 1 ),
                                   lepMinus->spParent( 1 ) );

        lac11 = EvtLeptonACurrent( lepPlus->spParent( 0 ),
                                   lepMinus->spParent( 0 ) );
        lac21 = EvtLeptonACurrent( lepPlus->spParent( 1 ),
                                   lepMinus->spParent( 0 ) );
        lac12 = EvtLeptonACurrent( lepPlus->spParent( 0 ),
                                   lepMinus->spParent( 1 ) );
        lac22 = EvtLeptonACurrent( lepPlus->spParent( 1 ),
                                   lepMinus->spParent( 1 ) );

        amp.vertex( 0, 0, CKM_factor * ( lvc11 * T1 + lac11 * T2 ) );
        amp.vertex( 0, 1, CKM_factor * ( lvc12 * T1 + lac12 * T2 ) );
        amp.vertex( 1, 0, CKM_factor * ( lvc21 * T1 + lac21 * T2 ) );
        amp.vertex( 1, 1, CKM_factor * ( lvc22 * T1 + lac22 * T2 ) );

    } else {
        if ( bbarmesons.contains( parentID ) ) {
            // The amplitude for the decay B -> K* ell^+ ell^-
            // (barb -> barq ell^+ ell^- transition)

            T1 = a_barb2barq * hatP + b_barb2barq * hatq;
            T2 = c * hatP + d * hatq;

            lvc11 = EvtLeptonVCurrent( lepPlus->spParent( 1 ),
                                       lepMinus->spParent( 1 ) );
            lvc21 = EvtLeptonVCurrent( lepPlus->spParent( 0 ),
                                       lepMinus->spParent( 1 ) );
            lvc12 = EvtLeptonVCurrent( lepPlus->spParent( 1 ),
                                       lepMinus->spParent( 0 ) );
            lvc22 = EvtLeptonVCurrent( lepPlus->spParent( 0 ),
                                       lepMinus->spParent( 0 ) );

            lac11 = EvtLeptonACurrent( lepPlus->spParent( 1 ),
                                       lepMinus->spParent( 1 ) );
            lac21 = EvtLeptonACurrent( lepPlus->spParent( 0 ),
                                       lepMinus->spParent( 1 ) );
            lac12 = EvtLeptonACurrent( lepPlus->spParent( 1 ),
                                       lepMinus->spParent( 0 ) );
            lac22 = EvtLeptonACurrent( lepPlus->spParent( 0 ),
                                       lepMinus->spParent( 0 ) );

            amp.vertex( 0, 0, conj( CKM_factor ) * ( lvc11 * T1 + lac11 * T2 ) );
            amp.vertex( 0, 1, conj( CKM_factor ) * ( lvc12 * T1 + lac12 * T2 ) );
            amp.vertex( 1, 0, conj( CKM_factor ) * ( lvc21 * T1 + lac21 * T2 ) );
            amp.vertex( 1, 1, conj( CKM_factor ) * ( lvc22 * T1 + lac22 * T2 ) );
        }

        else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "\n\n The function EvtbTosllScalarAmpNew::CalcAmp(...)"
                << "\n Wrong B-meson number" << std::endl;
            ::abort();
        }
    }
}

//
// The decays B -> P ell^+ ell^- maximum probability calculation for the
// d^2\Gamma/dq^2 d\cos\theta distribution.
//
// \theta - the angle between the final P-meson and ell^- directions in the
//          B-meson rest frame.
//
// If ias=0 (nonresonant case), the maximum is achieved at (s,t) plane!
// If ias=1 (resonant case), the maximum is achieved at q2=M^2_{J/\psi}.
//
double EvtbTosllScalarAmpNewExt::CalcMaxProb(
    EvtId parnum, EvtId mesnum, EvtId l1num, EvtId l2num,
    EvtbTosllFFNew* formFactors, EvtbTosllWilsCoeffNLO* WilsCoeff, double mu,
    int Nf, int res_swch, int ias, double CKM_A, double CKM_lambda,
    double CKM_barrho, double CKM_bareta, double ReA7, double ImA7,
    double ReA10, double ImA10 )
{
    double maxfoundprob = -100.0;    // maximum of the probability

    double M1 = EvtPDL::getMeanMass( parnum );    // B - meson mass
    double M2 = EvtPDL::getMeanMass( mesnum );    // P - meson mass
    double ml = EvtPDL::getMeanMass( l1num );     // leptonic mass

    if ( res_swch == 0 ) {
        double s, t_for_s;         // Mandelstam variables
        double s_min, s_max;       // s-variable boundaries
        double t_plus, t_minus;    // t-variable boundaries for current s-variable
        double ds, dt;
        int j, k;
        int max_j, max_k;

        s_min = 4.0 * pow( ml, 2.0 );       // minimum value of s-variable
        s_max = pow( ( M1 - M2 ), 2.0 );    // maximum value of s-variable

        max_j = 1000;
        ds = ( s_max - s_min ) / ( (double)max_j );
        if ( ds < 0.0 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "\n\n In the function EvtbTosllScalarAmpNew::CalcScalarMaxProb(...)"
                << "\n ds    = " << ds << " < 0."
                << "\n s_min = " << s_min << "\n s_max = " << s_max
                << "\n M1      = " << M1 << "\n M2      = " << M2
                << "\n ml      = " << ml << std::endl;
            ::abort();
        }

        // The maximum probability calculation
        // from s_min to s_max
        for ( j = max_j / 3; j < max_j; j++ ) {
            s = s_min + ds * ( (double)j );

            t_plus = pow( M1, 2.0 ) + pow( M2, 2.0 ) + 2.0 * pow( ml, 2.0 ) - s;
            t_plus = t_plus +
                     sqrt( 1.0 - 4.0 * pow( ml, 2.0 ) / s ) *
                         sqrt( lambda( s, pow( M1, 2.0 ), pow( M2, 2.0 ) ) );
            t_plus *= 0.5;

            t_minus = pow( M1, 2.0 ) + pow( M2, 2.0 ) + 2.0 * pow( ml, 2.0 ) - s;
            t_minus = t_minus -
                      sqrt( 1.0 - 4.0 * pow( ml, 2.0 ) / s ) *
                          sqrt( lambda( s, pow( M1, 2.0 ), pow( M2, 2.0 ) ) );
            t_minus *= 0.5;

            max_k = 1000;
            dt = ( t_plus - t_minus ) / ( (double)max_k );
            if ( fabs( dt ) < 0.00001 )
                dt = 0.0;

            if ( dt <= ( -0.00001 ) ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "\n\n In the function EvtbTosllScalarAmpNew::CalcScalarMaxProb(...)"
                    << "\n dt      = " << dt << " < 0."
                    << "\n s       = " << s << "\n s_min   = " << s_min
                    << "\n s_max   = " << s_max << "\n ds      = " << ds
                    << "\n j       = " << j << "\n t_plus  = " << t_plus
                    << "\n t_minus = " << t_minus << "\n M1      = " << M1
                    << "\n M2      = " << M2 << "\n ml      = " << ml
                    << std::endl;
                ::abort();
            }

            // from t_minus to t_plus
            for ( k = 0; k < max_k; k++ ) {
                t_for_s = t_minus + dt * ( (double)k );

                if ( ( t_for_s > t_plus ) && ( t_for_s <= ( 1.0001 * t_plus ) ) ) {
                    t_for_s = t_plus;
                }
                if ( t_for_s > ( 1.0001 * t_plus ) ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "\n\n In the function EvtbTosllScalarAmpNew::CalcScalarMaxProb(...)"
                        << "\n t_for_s = " << t_for_s
                        << " > t_plus = " << t_plus << " ! "
                        << "\n t_minus = " << t_minus << "\n dt      = " << dt
                        << "\n k       = " << k << "\n s       = " << s
                        << "\n M1      = " << M1 << "\n M2      = " << M2
                        << "\n ml      = " << ml << std::endl;
                    ::abort();
                }

                // B-meson rest frame particles and they kinematics inicialization
                double EV, El2;
                EV = ( pow( M1, 2.0 ) + pow( M2, 2.0 ) - s ) /
                     ( 2.0 * M1 );    // P-meson energy
                El2 = ( s + t_for_s - pow( M2, 2.0 ) - pow( ml, 2.0 ) ) /
                      ( 2.0 * M1 );    // ell^- energy

                double modV, modl2;
                modV = sqrt( pow( EV, 2.0 ) - pow( M2, 2.0 ) );
                modl2 = sqrt( pow( El2, 2.0 ) - pow( ml, 2.0 ) );

                double cosVellminus;    // angle between the P-meson and ell^- directions
                cosVellminus = ( pow( M2, 2.0 ) + pow( ml, 2.0 ) +
                                 2.0 * EV * El2 - t_for_s ) /
                               ( 2.0 * modV * modl2 );
                if ( ( fabs( cosVellminus ) > 1.0 ) &&
                     ( fabs( cosVellminus ) <= 1.0001 ) ) {
                    //            EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                    //               << "\n Debug in the function EvtbTosllScalarAmpNew::CalcMaxProb(...):"
                    //               << "\n cos(theta) = " << cosVellminus
                    //               << std::endl;
                    cosVellminus = cosVellminus / fabs( cosVellminus );
                }
                if ( ( modV <= 0.000001 ) || ( modl2 <= 0.000001 ) ) {
                    cosVellminus = cosVellminus / fabs( cosVellminus );
                    //            EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                    //               << "\n Debug in the function EvtbTosllScalarAmpNew::CalcMaxProb(...):"
                    //               << "\n modV       = " << modV
                    //               << "\n modl2      = " << modl2
                    //               << "\n cos(theta) = " << cosVellminus
                    //               << "\n s          = " << s
                    //               << "\n t_for_s    = " << t_for_s
                    //               << "\n s_min      = " << s_min
                    //               << "\n s_max      = " << s_max
                    //               << "\n t_plus     = " << t_plus
                    //               << "\n t_minus    = " << t_minus
                    //               << "\n dt         = " << dt
                    //               << "\n EV         = " << EV
                    //               << "\n El2        = " << El2
                    //               << "\n M2         = " << M2
                    //               << "\n ml         = " << ml
                    //               << std::endl;
                }
                if ( fabs( cosVellminus ) > 1.0001 ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "\n\n In the function EvtbTosllScalarAmpNew::CalcMaxProb(...)"
                        << "\n |cos(theta)| = " << fabs( cosVellminus ) << " > 1"
                        << "\n s       = " << s << "\n t_for_s = " << t_for_s
                        << "\n s_min = " << s_min << "\n s_max = " << s_max
                        << "\n t_plus  = " << t_plus << "\n t_minus = " << t_minus
                        << "\n dt      = " << dt << "\n EV      = " << EV
                        << "\n El2     = " << El2 << "\n modV    = " << modV
                        << "\n modl2   = " << modl2 << "\n M2      = " << M2
                        << "\n ml      = " << ml << std::endl;
                    ::abort();
                }

                EvtVector4R p1, p2, k1, k2;
                p1.set( M1, 0.0, 0.0, 0.0 );
                p2.set( EV, modV, 0.0, 0.0 );
                k2.set( El2, modl2 * cosVellminus,
                        -modl2 * sqrt( 1.0 - pow( cosVellminus, 2.0 ) ), 0.0 );
                k1 = p1 - p2 - k2;

                //          EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                //              << "\n Debug in the function EvtbTosllScalarAmpNew::CalcMaxProb(...):"
                //              << "\n mu =" << mu << " Nf =" << Nf
                //              << " res_swch =" << res_swch
                //              << " ias =" << ias
                //              << "\n M1 = " << M1
                //              << "\n M2 = " << M2
                //              << "\n ml = " << ml
                //              << "\n s  = " << s
                //              << "\n t_for_s = " << t_for_s
                //              << "\n EV      = " << EV
                //              << "\n El1     = " << El1
                //              << "\n El2     = " << El2
                //              << "\n modV    = " << modV
                //              << "\n modl1   = " << modl1
                //              << "\n modl2   = " << modl2
                //              << "\n cos(theta) = " << cosVellminus
                //              << "\n p1 =" << p1
                //              << "\n p2 =" << p2
                //              << "\n k1 =" << k1
                //              << "\n k2 =" << k2
                //              << std::endl;

                // B-meson state preparation at the rest frame of B-meson
                EvtScalarParticle* scalar_part;
                EvtParticle* root_part;
                scalar_part = new EvtScalarParticle;

                scalar_part->noLifeTime();
                scalar_part->init( parnum, p1 );
                root_part = (EvtParticle*)scalar_part;
                root_part->setDiagonalSpinDensity();

                // Amplitude initialization
                EvtId listdaug[3];
                listdaug[0] = mesnum;
                listdaug[1] = l1num;
                listdaug[2] = l2num;

                EvtAmp amp;
                amp.init( parnum, 3, listdaug );

                // Daughters states preparation at the rest frame of B-meson
                root_part->makeDaughters( 3, listdaug );

                EvtParticle *vect, *lep1, *lep2;
                vect = root_part->getDaug( 0 );
                lep1 = root_part->getDaug( 1 );
                lep2 = root_part->getDaug( 2 );

                vect->noLifeTime();
                lep1->noLifeTime();
                lep2->noLifeTime();

                vect->init( mesnum, p2 );
                lep1->init( l1num, k1 );
                lep2->init( l2num, k2 );

                EvtSpinDensity rho;
                rho.setDiag( root_part->getSpinStates() );

                // The amplitude calculation at the
                // "maximum amplitude" kinematical configuration
                CalcAmp( root_part, amp, formFactors, WilsCoeff, mu, Nf,
                         res_swch, ias, CKM_A, CKM_lambda, CKM_barrho,
                         CKM_bareta, ReA7, ImA7, ReA10, ImA10 );

                // Now find the probability at this q2 and cos theta lepton point
                double nikmax = rho.normalizedProb( amp.getSpinDensity() );

                if ( nikmax > maxfoundprob ) {
                    maxfoundprob = nikmax;
                    //             EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                    //                    << "\n maxfoundprob ( s =" << s << ",  t = " << t_for_s << " ) = "
                    //                    << maxfoundprob
                    //                    << "\n k =" << k
                    //                    << std::endl;
                }

                delete scalar_part;
                //          delete root_part;
                delete vect;
                delete lep1;
                delete lep2;

            }    // for(k=0; k<=max_k; k++)
        }        // for(j=0; j<=max_j; j++)

    }    // if(res_swch==0)

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if ( res_swch == 1 ) {
        double s, t_for_s;         // Mandelstam variables
        double t_plus, t_minus;    // t-variable boundaries for current s-variable
        double dt;
        int k;

        s = pow( 3.09688, 2.0 );    // s = (M_{J/\psi})^2

        t_plus = pow( M1, 2.0 ) + pow( M2, 2.0 ) + 2.0 * pow( ml, 2.0 ) - s;
        t_plus = t_plus + sqrt( 1.0 - 4.0 * pow( ml, 2.0 ) / s ) *
                              sqrt( lambda( s, pow( M1, 2.0 ), pow( M2, 2.0 ) ) );
        t_plus *= 0.5;

        t_minus = pow( M1, 2.0 ) + pow( M2, 2.0 ) + 2.0 * pow( ml, 2.0 ) - s;
        t_minus = t_minus -
                  sqrt( 1.0 - 4.0 * pow( ml, 2.0 ) / s ) *
                      sqrt( lambda( s, pow( M1, 2.0 ), pow( M2, 2.0 ) ) );
        t_minus *= 0.5;

        dt = ( t_plus - t_minus ) / 1000.0;

        // The maximum probability calculation
        for ( k = 0; k < 1000; k++ ) {
            t_for_s = t_minus + dt * ( (double)k );

            if ( ( t_for_s > t_plus ) && ( t_for_s <= ( 1.0001 * t_plus ) ) ) {
                t_for_s = t_plus;
            }
            if ( t_for_s > ( 1.0001 * t_plus ) ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "\n\n In the function EvtbTosllScalarAmpNew::CalcMaxProb(...)"
                    << "\n t_for_s = " << t_for_s << " > t_plus = " << t_plus
                    << " ! "
                    << "\n t_minus = " << t_minus << "\n dt      = " << dt
                    << "\n k       = " << k << "\n s       = " << s
                    << "\n M1      = " << M1 << "\n M2      = " << M2
                    << "\n ml      = " << ml << std::endl;
                ::abort();
            }

            // B-meson rest frame particles and they kinematics inicialization
            double EV, El2;
            EV = ( pow( M1, 2.0 ) + pow( M2, 2.0 ) - s ) /
                 ( 2.0 * M1 );    // V-meson energy
            El2 = ( s + t_for_s - pow( M2, 2.0 ) - pow( ml, 2.0 ) ) /
                  ( 2.0 * M1 );    // ell^- energy

            double modV, modl2;
            modV = sqrt( pow( EV, 2.0 ) - pow( M2, 2.0 ) );
            modl2 = sqrt( pow( El2, 2.0 ) - pow( ml, 2.0 ) );

            double cosVellminus;    // angle between the vector meson and ell^- directions
            cosVellminus = ( pow( M2, 2.0 ) + pow( ml, 2.0 ) + 2.0 * EV * El2 -
                             t_for_s ) /
                           ( 2.0 * modV * modl2 );
            if ( ( fabs( cosVellminus ) > 1.0 ) &&
                 ( fabs( cosVellminus ) <= 1.0001 ) ) {
                //         EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                //           << "\n Debug in the function EvtbTosllScalarAmpNew::CalcMaxProb(...):"
                //           << "\n cos(theta) = " << cosVellminus
                //           << std::endl;
                cosVellminus = cosVellminus / fabs( cosVellminus );
            }
            if ( ( modV <= 0.000001 ) || ( modl2 <= 0.000001 ) ) {
                cosVellminus = cosVellminus / fabs( cosVellminus );
                //         EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                //            << "\n Debug in the function EvtbTosllScalarAmpNew::CalcMaxProb(...):"
                //            << "\n modV       = " << modV
                //            << "\n modl2      = " << modl2
                //            << "\n cos(theta) = " << cosVellminus
                //            << "\n s          = " << s
                //            << "\n t_for_s    = " << t_for_s
                //            << "\n t_plus     = " << t_plus
                //            << "\n t_minus    = " << t_minus
                //            << "\n dt         = " << dt
                //            << "\n EV         = " << EV
                //            << "\n El2        = " << El2
                //            << "\n M2         = " << M2
                //            << "\n ml         = " << ml
                //            << std::endl;
            }
            if ( fabs( cosVellminus ) > 1.0001 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "\n\n In the function EvtbTosllScalarAmpNew::CalcMaxProb(...)"
                    << "\n |cos(theta)| = " << fabs( cosVellminus ) << " > 1"
                    << "\n s       = " << s << "\n t_for_s = " << t_for_s
                    << "\n t_plus  = " << t_plus << "\n t_minus = " << t_minus
                    << "\n dt      = " << dt << "\n EV      = " << EV
                    << "\n El2     = " << El2 << "\n modV    = " << modV
                    << "\n modl2   = " << modl2 << "\n M2      = " << M2
                    << "\n ml      = " << ml << std::endl;
                ::abort();
            }

            EvtVector4R p1, p2, k1, k2;
            p1.set( M1, 0.0, 0.0, 0.0 );
            p2.set( EV, modV, 0.0, 0.0 );
            k2.set( El2, modl2 * cosVellminus,
                    -modl2 * sqrt( 1.0 - pow( cosVellminus, 2.0 ) ), 0.0 );
            k1 = p1 - p2 - k2;

            //       EvtGenReport(EVTGEN_NOTICE,"EvtGen")
            //           << "\n Debug in the function EvtbTosllScalarAmpNew::CalcMaxProb(...):"
            //           << "\n mu =" << mu << " Nf =" << Nf
            //           << " res_swch =" << res_swch
            //           << " ias =" << ias
            //           << "\n M1 = " << M1
            //           << "\n M2 = " << M2
            //           << "\n ml = " << ml
            //           << "\n s  = " << s
            //           << "\n t_for_s = " << t_for_s
            //           << "\n EV      = " << EV
            //           << "\n El1     = " << El1
            //           << "\n El2     = " << El2
            //           << "\n modV    = " << modV
            //           << "\n modl1   = " << modl1
            //           << "\n modl2   = " << modl2
            //           << "\n cos(theta) = " << cosVellminus
            //           << "\n p1 =" << p1
            //           << "\n p2 =" << p2
            //           << "\n k1 =" << k1
            //           << "\n k2 =" << k2
            //           << std::endl;

            // B-meson state preparation at the rest frame of B-meson
            EvtScalarParticle* scalar_part;
            EvtParticle* root_part;
            scalar_part = new EvtScalarParticle;

            scalar_part->noLifeTime();
            scalar_part->init( parnum, p1 );
            root_part = (EvtParticle*)scalar_part;
            root_part->setDiagonalSpinDensity();

            // Amplitude initialization
            EvtId listdaug[3];
            listdaug[0] = mesnum;
            listdaug[1] = l1num;
            listdaug[2] = l2num;

            EvtAmp amp;
            amp.init( parnum, 3, listdaug );

            // Daughters states preparation at the rest frame of B-meson
            root_part->makeDaughters( 3, listdaug );

            EvtParticle *vect, *lep1, *lep2;
            vect = root_part->getDaug( 0 );
            lep1 = root_part->getDaug( 1 );
            lep2 = root_part->getDaug( 2 );

            vect->noLifeTime();
            lep1->noLifeTime();
            lep2->noLifeTime();

            vect->init( mesnum, p2 );
            lep1->init( l1num, k1 );
            lep2->init( l2num, k2 );

            EvtSpinDensity rho;
            rho.setDiag( root_part->getSpinStates() );

            // The amplitude calculation at the
            // "maximum amplitude" kinematical configuration
            CalcAmp( root_part, amp, formFactors, WilsCoeff, mu, Nf, res_swch,
                     ias, CKM_A, CKM_lambda, CKM_barrho, CKM_bareta, ReA7, ImA7,
                     ReA10, ImA10 );

            // Now find the probability at this q2 and cos theta lepton point
            double nikmax = rho.normalizedProb( amp.getSpinDensity() );

            if ( nikmax > maxfoundprob ) {
                maxfoundprob = nikmax;
                //          EvtGenReport(EVTGEN_NOTICE,"EvtGen")
                //                 << "\n maxfoundprob ( s =" << s << ",  t = " << t_for_s << " ) = "
                //                 << maxfoundprob
                //                 << "\n k =" << k
                //                 << std::endl;
            }

            delete scalar_part;
            //       delete root_part;
            delete vect;
            delete lep1;
            delete lep2;

        }    // for(k=0; k<=1000; k++)

    }    // if(res_swch==1)

    if ( maxfoundprob <= 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n In the function EvtbTosllScalarAmpNew::CalcMaxProb(...)"
            << "\n maxfoundprob = " << maxfoundprob << " <0 or =0!"
            << "\n res_swch     = " << res_swch << std::endl;
        ::abort();
    }

    EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
        << "\n maxfoundprob (...) = " << maxfoundprob << std::endl;

    maxfoundprob *= 1.01;

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen")
    //         << "\n ***************************************************************************"
    //         << "\n The function EvtbTosllScalarAmpNew::CalcMaxProb(...) passed with arguments:"
    //         << "\n mu =" << mu << " Nf =" << Nf
    //         << " res_swch =" << res_swch
    //         << " ias =" << ias
    //         << " \n s_at_max = " << s_at_max
    //         << " t_at_max = " << t_at_max
    //         << "\n The distribution maximum maxfoundprob =" << maxfoundprob
    //         << "\n ***************************************************************************"
    //         << std::endl;

    return maxfoundprob;
}

// Triangular function
double EvtbTosllScalarAmpNewExt::lambda( double a, double b, double c )
{
    double l;

    l = pow( a, 2.0 ) + pow( b, 2.0 ) + pow( c, 2.0 ) - 2.0 * a * b -
        2.0 * a * c - 2.0 * b * c;

    return l;
}
