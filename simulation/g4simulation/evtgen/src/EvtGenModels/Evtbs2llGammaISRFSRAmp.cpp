
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

#include "EvtGenModels/Evtbs2llGammaISRFSRAmp.hh"

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
#include "EvtGenBase/EvtVectorParticle.hh"

#include "EvtGenModels/EvtbTosllWilsCoeffNLO.hh"
#include "EvtGenModels/Evtbs2llGammaFFMNT.hh"

#include <cstdlib>

// input:   *parent      - the pointer to the parent particle (B-meson, the
//                                          object of the EvtParticle class);
//          *formFactors - the pointer to instance of EvtbTosllGammaFF class object;
//          *WilsCoeff   - the pointer to the Standart Model Wilson Coefficients class;
//           mu          - the scale parameter, GeV;
//           Nf          - number of "effective" flavors (for b-quark Nf=5);
//           sr          - state radiation type:
//                         = 0 the ISR only,
//                         = 1 the FSR only,
//                         = 2 both ISR + FSR (== BSTOGLLMNT model);
//           res_swch    - resonant switching parameter:
//                         = 0 the resonant contribution switched OFF,
//                         = 1 the resonant contribution switched ON;
//           ias         - switching parameter for \alpha_s(M_Z) value:
//                         = 0 PDG 1sigma minimal alpha_s(M_Z),
//                         = 1 PDG average value  alpha_s(M_Z),
//                         = 2 PDG 1sigma maximal alpha_s(M_Z).
//           Egamma_min  - photon energy cut, GeV;
//           Wolfenstein parameterization for CKM matrix
//                         CKM_A, CKM_lambda, CKM_barrho, CKM_bareta

void Evtbs2llGammaISRFSRAmp::CalcAmp( EvtParticle* parent, EvtAmp& amp,
                                      Evtbs2llGammaFF* formFactors,
                                      EvtbTosllWilsCoeffNLO* WilsCoeff,
                                      double mu, int Nf, int sr, int res_swch,
                                      int ias, double Egamma_min, double CKM_A,
                                      double CKM_lambda, double CKM_barrho,
                                      double CKM_bareta, double mumumass_min )
{
    //  FILE *mytest;

    int iG = 0;    // photon is the first daughter particle
    int il1 = 1,
        il2 = 2;    // leptons are the second and thirds daughter particles

    EvtComplex unit1( 1.0, 0.0 );    // real unit
    EvtComplex uniti( 0.0, 1.0 );    // imaginary unit

    double M1 = parent->mass();                    // B - meson mass, GeV
    double ml = parent->getDaug( il1 )->mass();    // leptonic mass, GeV
    double mq = 0.0;    // light quark mass from the dispersion QM, GeV
    double mc = formFactors->getQuarkMass(
        4 );    // m_c mass from the dispersion QM, GeV
    double mb = formFactors->getQuarkMass(
        5 );               // m_b mass from the dispersion QM, GeV
    double Mw = 80.403;    // GeV W-boson mass, GeV
    double mt = 174.2;     // GeV t-quark mass, GeV
    double fb = 0.0;       // leptonic decay constant of B-meson, Gev

    EvtComplex Vtb, Vtq, Vub, Vuq, Vcb,
        Vcq;                  // V_{tb}, V_{tq}, V_{ub}, V_{uq}, V_{cb}, V_{cq}
    EvtComplex CKM_factor;    // V^*_{tq}*V_{tb}, where q={d,s}
    EvtComplex lambda_qu;     // V^*_{uq}*V_{ub}/V^*_{tq}*V_{tb}, where q={d,s}
    EvtComplex lambda_qc;     // V^*_{cq}*V_{cb}/V^*_{tq}*V_{tb}, where q={d,s}
    double Relambda_qu, Imlambda_qu;

    // to find charges of ell^+ and ell^- in the B-meson  daughters
    int charge1 = EvtPDL::chg3( parent->getDaug( il1 )->getId() );
    int charge2 = EvtPDL::chg3( parent->getDaug( il2 )->getId() );
    if ( charge1 == 0 || charge2 == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n The function EvtbsTollGammaAmp::CalcAmp(...)"
            << "\n Error in the leptonic charge getting!"
            << "\n charge1             =" << charge1
            << "\n charge2             =" << charge2 << "\n charge gamma        ="
            << EvtPDL::chg3( parent->getDaug( iG )->getId() )
            << "\n number of daughters =" << parent->getNDaug() << std::endl;
        ::abort();
    }

    EvtParticle* lepPlus = 0;
    EvtParticle* lepMinus = 0;

    lepPlus = ( charge1 > charge2 )
                  ? parent->getDaug( il1 )
                  : parent->getDaug( il2 );    // positive charged
    lepMinus = ( charge1 < charge2 )
                   ? parent->getDaug( il1 )
                   : parent->getDaug( il2 );    // negative charged

    EvtVector4R p = parent->getP4Restframe();    // B-meson momentum in the B-rest frame
    EvtVector4R k =
        parent->getDaug( iG )->getP4();    // 4-momentum of photon in the B-rest frame
    EvtVector4R q = p - k;    // transition 4-momentum q=p-k in the B-rest frame
    EvtVector4R p_1;          // 4-momentum of ell^+ in the B-rest frame
    EvtVector4R p_2;          // 4-momentum of ell^- in the B-rest frame

    // the preparation of the leptonic 4-momentums in the B-rest frame
    if ( charge1 > charge2 ) {
        p_1 = parent->getDaug( il1 )->getP4();
        p_2 = parent->getDaug( il2 )->getP4();
    } else {
        p_1 = parent->getDaug( il2 )->getP4();
        p_2 = parent->getDaug( il1 )->getP4();
    }

    EvtVector4R p_minus_p_1 =
        p - p_1;    // transition momentum of the B-meson and antilepton p-p_1
    EvtVector4R p_minus_p_2 =
        p - p_2;    // transition momentum of the B-meson and lepton p-p_2

    double q2 = q.mass2();             // Mandelstam variable s=q^2
    double p2 = p.mass2();             // p^2=M1^2
    double t = p_minus_p_1.mass2();    // Mandelstam variable t=(p-p_1)^2
    double u = p_minus_p_2.mass2();    // Mandelstam variable u=(p-p_2)^2

    // scalar products
    double pk = 0.5 * ( p2 - q2 );                // (p*k)
    double p1k = 0.5 * ( pow( ml, 2.0 ) - u );    // (p1*k)
    double p2k = 0.5 * ( pow( ml, 2.0 ) - t );    // (p2*k)

    double hatq2 = q2 / ( M1 * M1 );    // \hat s = q^2/M_1^2
    double Egam = 0.5 * M1 *
                  ( 1 - hatq2 );    // photon energy in the B-meson rest frame

    EvtVector4R hatp = p / M1;
    EvtVector4R hatk = k / M1;

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen")
    //         << "\n\n The function EvtbsTollGammaISRFSRAmp::CalcAmp(...)"
    //         << "\n q =  p-k  =" << p-k     << "   q^2 = " << (p-k).mass2()
    //         << "\n q = p1+p2 =" << p_1+p_2 << "   q^2 = " << (p_1+p_2).mass2()
    //         << "\n m_ell     =" << parent->getDaug(il1)->mass()
    //         << "\n m_ell     =" << parent->getDaug(il2)->mass()
    //         << "\n m_gamma   =" << parent->getDaug(iG)->mass()
    //         << std::endl;

    EvtId idparent = parent->getId();    // B-meson Id

    if ( idparent == EvtPDL::getId( std::string( "B_s0" ) ) ||
         idparent == EvtPDL::getId( std::string( "anti-B_s0" ) ) ) {
        mq = formFactors->getQuarkMass( 3 );    // m_s mass from the dispersion QM
        fb = 0.24;                              // leptonic decay constant
            // V_{ts}
        Vtq = unit1 * ( 1.0 - 0.5 * pow( CKM_lambda, 2.0 ) ) +
              pow( CKM_lambda, 2.0 ) *
                  ( CKM_barrho * unit1 + CKM_bareta * uniti ) /
                  sqrt( 1.0 - pow( CKM_lambda, 2.0 ) );
        Vtq = -CKM_A * pow( CKM_lambda, 2.0 ) * Vtq;
        // V_{us}
        Vuq = CKM_lambda * unit1;
        // V_{cs}
        Vcq = unit1 - 0.5 * pow( CKM_lambda, 2.0 ) -
              0.125 * pow( CKM_lambda, 4.0 ) * ( 1.0 + 4.0 * pow( CKM_A, 2.0 ) );
    }

    if ( idparent == EvtPDL::getId( std::string( "B0" ) ) ||
         idparent == EvtPDL::getId( std::string( "anti-B0" ) ) ) {
        mq = formFactors->getQuarkMass( 2 );    // m_d mass from the dispersion QM
        fb = 0.20;                              // leptonic decay constant
            // V_{td}
        Vtq = unit1 - ( 1.0 - 0.5 * pow( CKM_lambda, 2.0 ) ) *
                          ( CKM_barrho * unit1 + CKM_bareta * uniti ) /
                          sqrt( 1.0 - pow( CKM_lambda, 2.0 ) );
        Vtq = CKM_A * pow( CKM_lambda, 3.0 ) * Vtq;
        // V_{ud}
        Vuq = unit1 * ( 1.0 - 0.5 * pow( CKM_lambda, 2.0 ) -
                        0.125 * pow( CKM_lambda, 4.0 ) );
        // V_{cd}
        Vcq = unit1 *
              ( -CKM_lambda +
                0.5 * pow( CKM_A, 2.0 ) * pow( CKM_lambda, 5.0 ) *
                    ( 1.0 - 2.0 * ( CKM_barrho * unit1 + CKM_bareta * uniti ) /
                                sqrt( 1.0 - pow( CKM_lambda, 2.0 ) ) ) );
    }

    if ( mq < 0.001 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n The function EvtbsTollGammaISRFSRAmp::CalcAmp(..// 4-momentum of ell^+.)"
            << "\n Error in the model set!"
            << " mq = " << mq << std::endl;
        ::abort();
    }

    Vtb = unit1 * ( 1.0 - 0.5 * pow( CKM_A * CKM_lambda * CKM_lambda,
                                     2.0 ) );    // V_{tb}
    Vub = CKM_A * pow( CKM_lambda, 3.0 ) *
          ( CKM_barrho * unit1 - CKM_bareta * uniti ) /
          sqrt( 1.0 - pow( CKM_lambda, 2.0 ) );      // V_{ub}
    Vcb = unit1 * CKM_A * pow( CKM_lambda, 2.0 );    // V_{cb}

    CKM_factor = conj( Vtq ) * Vtb;    // V^*_{tq}*V_{tb}

    lambda_qu = conj( Vuq ) * Vub /
                CKM_factor;    // V^*_{uq}*V_{ub}/V^*_{tq}*V_{tb}
    Relambda_qu = real( lambda_qu );
    Imlambda_qu = imag( lambda_qu );

    lambda_qc = conj( Vcq ) * Vcb /
                CKM_factor;    // V^*_{cq}*V_{cb}/V^*_{tq}*V_{tb}

    // The Wilson Coefficients preparation according to the paper
    // A.J.Buras, M.Munz, Phys.Rev.D52, p.189 (1995)
    double c1, c2;
    EvtComplex a1, c7gam, c9eff_b2q, c9eff_barb2barq, c10a;

    // foton energy cut and removal of the J/psi amd psi' resonant area
    if ( Egam < Egamma_min || ( res_swch == 1 && q2 >= 9.199 && q2 <= 15.333 ) ||
         ( q2 <= mumumass_min * mumumass_min ) ) {
        c1 = 0.0;
        c2 = 0.0;
        a1 = unit1 * 0.0;
        c7gam = unit1 * 0.0;
        c9eff_b2q = unit1 * 0.0;
        c9eff_barb2barq = unit1 * 0.0;
        c10a = unit1 * 0.0;
    } else {
        c1 = WilsCoeff->C1( mu, Mw, Nf, ias );
        c2 = WilsCoeff->C2( mu, Mw, Nf, ias );
        a1 = unit1 * ( c1 + c2 / 3.0 );
        c7gam = WilsCoeff->GetC7Eff( mu, Mw, mt, Nf, ias );
        c9eff_b2q = WilsCoeff->GetC9Eff( 0, res_swch, ias, Nf, q2, mb, mq, mc, mu,
                                         mt, Mw, ml, Relambda_qu, Imlambda_qu );
        c9eff_barb2barq = WilsCoeff->GetC9Eff( 1, res_swch, ias, Nf, q2, mb, mq,
                                               mc, mu, mt, Mw, ml, Relambda_qu,
                                               Imlambda_qu );
        c10a = WilsCoeff->GetC10Eff( mt, Mw );
    }

    EvtComplex Fv,
        Fa;    // The change of the sign is included in the amplitudes definition!
    EvtComplex Ftv_b2q, Ftv_barb2barq;
    EvtComplex Fta_b2q, Fta_barb2barq;

    // foton energy cut and removal of the J/psi amd psi' resonant area
    if ( Egam < Egamma_min || ( res_swch == 1 && q2 >= 9.199 && q2 <= 15.333 ) ||
         ( q2 <= mumumass_min * mumumass_min ) ) {
        fb = 0.0;
        Fa = unit1 * 0.0;
        Fv = unit1 * 0.0;
        Fta_b2q = unit1 * 0.0;
        Fta_barb2barq = unit1 * 0.0;
        Ftv_b2q = unit1 * 0.0;
        Ftv_barb2barq = unit1 * 0.0;
    } else {
        if ( fb < 0.01 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "\n\n The function EvtbsTollGammaISRFSRAmp::CalcAmp(...)"
                << "\n Leptonic decay constant fb is not uninitialized in this function!"
                << " fb = " << fb << std::endl;
            ::abort();
        }

        // For \bar B^0_q -> l^+ l^- gamma
        formFactors->getPhotonFF( 0, fb, parent->getId(), q2, M1, mb, mq, c7gam,
                                  a1, lambda_qu, lambda_qc, Fv, Fa, Ftv_b2q,
                                  Fta_b2q );

        // For B^0_q -> l^+ l^- gamma
        formFactors->getPhotonFF( 1, fb, parent->getId(), q2, M1, mb, mq, c7gam,
                                  a1, lambda_qu, lambda_qc, Fv, Fa,
                                  Ftv_barb2barq, Fta_barb2barq );
    }

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen")
    //      << "\n ============================================================================"
    //      << "\n ============================================================================"
    //      << "\n\n The function Evtbs2llGammaISRFSRAmp::CalcAmp(...) passed."
    //      << "\n Particle masses:"
    //      << "\n B - meson mass M1 = " << M1
    //      << "\n photon minimum E  = " << Egamma_min
    //      << "\n q2                = " << q2
    //      << "\n leptonic mass  ml = " << ml
    //      << "\n light quark mass  = " << mq
    //      << "\n c - quark mass mc = " << mc
    //      << "\n b - quark mass mb = " << mb
    //      << "\n t - quark mass mt = " << mt
    //      << "\n W - boson mass Mw = " << Mw
    //      << "\n ============================================================================"
    //      << "\n Input parameters:"
    //      << "\n scale parameter         mu = " << mu
    //      << "\n number of flavors       Nf = " << Nf
    //      << "\n state radiation switching  = " << sr
    //      << "\n resonant switching         = " << res_swch
    //      << "\n parameter for alpha_s(M_Z) = " << ias
    //      << "\n photon energy cut (GeV)    = " << Egamma_min
    //      << "\n ============================================================================"
    //      << "\n Form-factors"
    //      << "\n Egam          = " << Egam
    //      << "\n Egamma_min    = " << Egamma_min
    //      << "\n Fv            = " << Fv
    //      << "\n Fa            = " << Fa
    //      << "\n Ftv_b2q       = " << Ftv_b2q
    //      << "\n Fta_b2q       = " << Fta_b2q
    //      << "\n Ftv_barb2barq = " << Ftv_barb2barq
    //      << "\n Fta_barb2barq = " << Fta_barb2barq
    //      << "\n ============================================================================"
    //      << "\n Wilson Coefficients:"
    //      << "\n Egam                = " << Egam
    //      << "\n Egamma_min          = " << Egamma_min
    //      << "\n Re(c7gam)           = " << real(c7gam)
    //        << " Im(c7gam)           = " << imag(c7gam)
    //      << "\n Re(c9eff_b2q)       = " << real(c9eff_b2q)
    //        << " Im(c9eff_b2q)       = " << imag(c9eff_b2q)
    //      << "\n Re(c9eff_barb2barq) = " << real(c9eff_barb2barq)
    //        << " Im(c9eff_barb2barq) = " << imag(c9eff_barb2barq)
    //      << "\n Re(c10a)            = " << real(c10a)
    //        << " Im(c10a)            = " << imag(c10a)
    //      << std::endl;

    // Hadronic matrix element coefficients
    EvtComplex a_b2q, a_barb2barq, b_b2q, b_barb2barq, e_b2q, e_barb2barq,
        f_b2q, f_barb2barq;
    EvtComplex brammS, brammT;

    a_b2q = c9eff_b2q * Fv + 2.0 * c7gam * Ftv_b2q * mb * M1 / q2;
    a_barb2barq = c9eff_barb2barq * Fv +
                  2.0 * c7gam * Ftv_barb2barq * mb * M1 / q2;

    b_b2q = ( c9eff_b2q * Fa + 2.0 * c7gam * Fta_b2q * mb * M1 / q2 ) * pk /
            ( M1 * M1 );
    b_barb2barq = ( c9eff_barb2barq * Fa +
                    2.0 * c7gam * Fta_barb2barq * mb * M1 / q2 ) *
                  pk / ( M1 * M1 );

    e_b2q = c10a * Fv;
    e_barb2barq = e_b2q;

    f_b2q = c10a * Fa * pk / ( M1 * M1 );
    f_barb2barq = f_b2q;

    brammS = 0.0;    // in the Bq-meson rest frame!
    brammT = 0.5 * c10a * ml * fb *
             ( 1.0 / p2k + 1.0 / p1k );    // for Bramsstrahlung

    // The separation of the ISR and FSR contributions
    if ( sr == 0 ) {    // ISR only
        brammS = 0.0;
        brammT = 0.0;
    }
    if ( sr == 1 ) {    // FSR only
        a_b2q = 0.0;
        a_barb2barq = 0.0;
        b_b2q = 0.0;
        b_barb2barq = 0.0;
        e_b2q = 0.0;
        e_barb2barq = 0.0;
        f_b2q = 0.0;
        f_barb2barq = 0.0;
    }

    EvtTensor4C T1, T2;    // hadronic matrix element tensor structures
    EvtVector4C E1, E2;
    EvtComplex E3;

    EvtVector4C epsG;    // photon polarisation vector

    int i;    // photon polarisations counter

    EvtVector4C lvc11, lvc12;    // spin structures for
    EvtVector4C lvc21, lvc22;    // the leptonic vector current

    EvtVector4C lac11, lac12;    // spin structures for
    EvtVector4C lac21, lac22;    // the leptonic axial current

    EvtComplex lsc11, lsc12;    // spin structures for
    EvtComplex lsc21, lsc22;    // the leptonic scalar current

    EvtTensor4C ltc11, ltc12;    // spin structures for
    EvtTensor4C ltc21, ltc22;    // the leptonic tensor current

    // B - and barB - mesons descriptors
    static EvtIdSet bmesons( "anti-B0", "anti-B_s0" );
    static EvtIdSet bbarmesons( "B0", "B_s0" );

    EvtId parentID = parent->getId();

    if ( bmesons.contains( parentID ) ) {
        // The amplitude for the decay barB -> gamma ell^+ ell^-  or
        // b \bar q -> gamma ell^+ ell^-

        T1 = -a_b2q * unit1 * dual( EvtGenFunctions::directProd( hatp, hatk ) ) -
             b_b2q * uniti * EvtTensor4C::g();

        T2 = -e_b2q * unit1 * dual( EvtGenFunctions::directProd( hatp, hatk ) ) -
             f_b2q * uniti * EvtTensor4C::g();

        // spin combinations for vector lepton current
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

        lsc11 = EvtLeptonSCurrent( lepPlus->spParent( 0 ),
                                   lepMinus->spParent( 0 ) );
        lsc21 = EvtLeptonSCurrent( lepPlus->spParent( 1 ),
                                   lepMinus->spParent( 0 ) );
        lsc12 = EvtLeptonSCurrent( lepPlus->spParent( 0 ),
                                   lepMinus->spParent( 1 ) );
        lsc22 = EvtLeptonSCurrent( lepPlus->spParent( 1 ),
                                   lepMinus->spParent( 1 ) );

        // \epsilon^{\alpha\beta\mu\nu}*TCurrent_{\mu\nu}
        ltc11 = dual( EvtLeptonTCurrent( lepPlus->spParent( 0 ),
                                         lepMinus->spParent( 0 ) ) );
        ltc21 = dual( EvtLeptonTCurrent( lepPlus->spParent( 1 ),
                                         lepMinus->spParent( 0 ) ) );
        ltc12 = dual( EvtLeptonTCurrent( lepPlus->spParent( 0 ),
                                         lepMinus->spParent( 1 ) ) );
        ltc22 = dual( EvtLeptonTCurrent( lepPlus->spParent( 1 ),
                                         lepMinus->spParent( 1 ) ) );

        // summing up photon polarisations
        for ( i = 0; i < 2; i++ ) {
            // conjaction of epsG (photon polarization vector)
            EvtVector4C epsG = parent->getDaug( 0 )->epsParentPhoton( i ).conj();

            // de-escalation T with epsG
            E1 = T1.cont2( epsG );
            E2 = T2.cont2( epsG );
            E3 = ( epsG * hatp ) * brammS;

            // foton energy cut and removal of the J/psi amd psi' resonant area
            if ( Egam < Egamma_min ||
                 ( res_swch == 1 && q2 >= 9.199 && q2 <= 15.333 ) ||
                 ( q2 <= mumumass_min * mumumass_min ) ) {
                CKM_factor = 0.0 * unit1;
            }

            // 1
            amp.vertex( i, 0, 0,
                        CKM_factor *
                            ( lvc11 * E1 + lac11 * E2 + uniti * lsc11 * E3 +
                              uniti * ( ( ltc11.cont2( hatp ) ) * epsG ) *
                                  brammT ) );
            //      EvtGenReport(EVTGEN_NOTICE,"EvtGen")
            //        << "\n 1" << CKM_factor*(lvc11*E1+lac11*E2+uniti*lsc11*E3+uniti*((ltc11.cont2(hatp))*epsG)*brammT)
            //        << std::endl;

            // 2
            amp.vertex( i, 0, 1,
                        CKM_factor *
                            ( lvc12 * E1 + lac12 * E2 + uniti * lsc12 * E3 +
                              uniti * ( ( ltc12.cont2( hatp ) ) * epsG ) *
                                  brammT ) );
            //      EvtGenReport(EVTGEN_NOTICE,"EvtGen")
            //        << "\n 2" << CKM_factor*(lvc12*E1+lac12*E2+uniti*lsc12*E3+uniti*((ltc12.cont2(hatp))*epsG)*brammT)
            //        << std::endl;

            // 3
            amp.vertex( i, 1, 0,
                        CKM_factor *
                            ( lvc21 * E1 + lac21 * E2 + uniti * lsc21 * E3 +
                              uniti * ( ( ltc21.cont2( hatp ) ) * epsG ) *
                                  brammT ) );
            //      EvtGenReport(EVTGEN_NOTICE,"EvtGen")
            //        << "\n 3" << CKM_factor*(lvc21*E1+lac21*E2+uniti*lsc21*E3+uniti*((ltc21.cont2(hatp))*epsG)*brammT)
            //        << std::endl;

            // 4
            amp.vertex( i, 1, 1,
                        CKM_factor *
                            ( lvc22 * E1 + lac22 * E2 + uniti * lsc22 * E3 +
                              uniti * ( ( ltc22.cont2( hatp ) ) * epsG ) *
                                  brammT ) );
            //      EvtGenReport(EVTGEN_NOTICE,"EvtGen")
            //        << "\n 4" << CKM_factor*(lvc22*E1+lac22*E2+uniti*lsc22*E3+uniti*((ltc22.cont2(hatp))*epsG)*brammT)
            //        << std::endl;
        }

    } else {
        if ( bbarmesons.contains( parentID ) ) {
            // The amplitude for the decay B -> gamma ell^+ ell^- or
            // q bar b -> gamma ell^+ ell^-

            T1 = -a_barb2barq * unit1 *
                     dual( EvtGenFunctions::directProd( hatp, hatk ) ) +
                 b_barb2barq * uniti * EvtTensor4C::g();

            T2 = -e_barb2barq * unit1 *
                     dual( EvtGenFunctions::directProd( hatp, hatk ) ) +
                 f_barb2barq * uniti * EvtTensor4C::g();

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

            lsc11 = EvtLeptonSCurrent( lepPlus->spParent( 1 ),
                                       lepMinus->spParent( 1 ) );
            lsc21 = EvtLeptonSCurrent( lepPlus->spParent( 0 ),
                                       lepMinus->spParent( 1 ) );
            lsc12 = EvtLeptonSCurrent( lepPlus->spParent( 1 ),
                                       lepMinus->spParent( 0 ) );
            lsc22 = EvtLeptonSCurrent( lepPlus->spParent( 0 ),
                                       lepMinus->spParent( 0 ) );

            // \epsilon^{\alpha\beta\mu\nu}*TCurrent_{\mu\nu}
            ltc11 = dual( EvtLeptonTCurrent( lepPlus->spParent( 1 ),
                                             lepMinus->spParent( 1 ) ) );
            ltc21 = dual( EvtLeptonTCurrent( lepPlus->spParent( 0 ),
                                             lepMinus->spParent( 1 ) ) );
            ltc12 = dual( EvtLeptonTCurrent( lepPlus->spParent( 1 ),
                                             lepMinus->spParent( 0 ) ) );
            ltc22 = dual( EvtLeptonTCurrent( lepPlus->spParent( 0 ),
                                             lepMinus->spParent( 0 ) ) );

            // summing up photon polarisations
            for ( i = 0; i < 2; i++ ) {
                EvtVector4C barepsG = parent->getDaug( 0 )->epsParentPhoton( i );

                E1 = T1.cont2( barepsG );
                E2 = T2.cont2( barepsG );
                E3 = ( barepsG * hatp ) * brammS;

                // foton energy cut and removal of the J/psi amd psi' resonant area
                if ( Egam < Egamma_min ||
                     ( res_swch == 1 && q2 >= 9.199 && q2 <= 15.333 ) ||
                     ( q2 <= mumumass_min * mumumass_min ) ) {
                    CKM_factor = 0.0 * unit1;
                }

                amp.vertex(
                    i, 1, 1,
                    conj( CKM_factor ) *
                        ( lvc11 * E1 + lac11 * E2 + uniti * lsc11 * E3 +    // -?
                          uniti * ( ( ltc11.cont2( hatp ) ) * epsG ) * brammT ) );
                amp.vertex(
                    i, 1, 0,
                    conj( CKM_factor ) *
                        ( lvc12 * E1 + lac12 * E2 + uniti * lsc12 * E3 +    // -?
                          uniti * ( ( ltc12.cont2( hatp ) ) * epsG ) * brammT ) );
                amp.vertex(
                    i, 0, 1,
                    conj( CKM_factor ) *
                        ( lvc21 * E1 + lac21 * E2 + uniti * lsc21 * E3 +    // -?
                          uniti * ( ( ltc21.cont2( hatp ) ) * epsG ) * brammT ) );
                amp.vertex(
                    i, 0, 0,
                    conj( CKM_factor ) *
                        ( lvc22 * E1 + lac22 * E2 + uniti * lsc22 * E3 +    // -?
                          uniti * ( ( ltc22.cont2( hatp ) ) * epsG ) * brammT ) );
            }

        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "\n\n The function Evtbs2llGammaISRFSRAmp::CalcAmp(...)"
                << "\n Wrong B-meson number" << std::endl;
            ::abort();
        }
    }
}

//
// The decays B -> Gamma ell^+ ell^- maximum probability calculation for the
// d^2\Gamma/dq^2 d\cos\theta distribution.
//
// \theta - the angle between the photon and ell^- directions in the
//          B-meson rest frame.
//
// If ias=0 (nonresonant case), the maximum is achieved at q2
//           B0s:                q2 = 4*ml^2, Mphi^2, q_max^2;
//           B0d:                q2 = 4*ml^2, Mrho^2, Momega^2, q_max^2;
// If ias=1 (resonant case), the maximum in the same points, because the
//          resonat area is remove
//
double Evtbs2llGammaISRFSRAmp::CalcMaxProb(
    EvtId parnum, EvtId photnum, EvtId l1num, EvtId l2num,
    Evtbs2llGammaFF* formFactors, EvtbTosllWilsCoeffNLO* WilsCoeff, double mu,
    int Nf, int sr, int res_swch, int ias, double Egamma_min, double CKM_A,
    double CKM_lambda, double CKM_barrho, double CKM_bareta, double mumumass_min )
{
    double maxfoundprob = -100.0;    // maximum of the probability

    double M1 = EvtPDL::getMeanMass( parnum );    // B - meson mass
    double ml = EvtPDL::getMeanMass( l1num );     // leptonic mass

    double Mrho = EvtPDL::getMeanMass(
        EvtPDL::getId( std::string( "rho0" ) ) );    // mass of the rho-meson, MeV
    double Momega = EvtPDL::getMeanMass( EvtPDL::getId(
        std::string( "omega" ) ) );    // mass of the omega-meson, MeV
    double Mphi = EvtPDL::getMeanMass(
        EvtPDL::getId( std::string( "phi" ) ) );    // mass of the phi-meson, MeV

    // EvtGenReport(EVTGEN_NOTICE,"EvtGen")
    //   << "\n M1         = " << M1
    //   << "\n ml         = " << ml
    //   << "\n Mrho       = " << Mrho
    //   << "\n Momega     = " << Momega
    //   << "\n Mphi       = " << Mphi
    //   << "\n Egamma_min = " << Egamma_min
    //   << std::endl;

    double list_of_max_q2_points[5];
    list_of_max_q2_points[0] = pow( 2.0 * ml, 2.0 );
    list_of_max_q2_points[1] = pow( Mrho, 2.0 );
    list_of_max_q2_points[2] = pow( Momega, 2.0 );
    list_of_max_q2_points[3] = pow( Mphi, 2.0 );
    list_of_max_q2_points[4] =
        pow( M1, 2.0 ) - 2.0 * M1 * Egamma_min;    // q^2_max at photon energy cut

    //  if(list_of_max_points[4]<0){
    //     EvtGenReport(EVTGEN_ERROR,"EvtGen")
    //       << "\n\n In the function EvtbsTollGammaAmp::CalcScalarMaxProb(...)"
    //       << "\n Bad photon energy cut: Egamma_min > M1 in the rest frame of B-meson!"
    //       << "\n q2_max     = " << list_of_max_points[4]
    //       << "\n M1         = " << M1
    //       << "\n Egamma_min = " << Egamma_min
    //       << std::endl;
    //  ::abort();
    //  }

    if ( Egamma_min > Mrho ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n In the function Evtbs2llGammaISRFSRAmp::CalcMaxProb(...)"
            << "\n Bad photon energy cut: Egamma_min > M_rho0 in the rest frame of B-meson."
            << "\n Mrho       = " << Mrho << "\n Egamma_min = " << Egamma_min
            << std::endl;
        ::abort();
    }

    if ( Egamma_min <= 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n In the function Evtbs2llGammaISRFSRAmp::CalcMaxProb(...)"
            << "\n Bad photon energy cut: Egamma_min <= 0 in the rest frame of B-meson."
            << "\n Egamma_min = " << Egamma_min << std::endl;
        ::abort();
    }

    if ( res_swch == 0 || res_swch == 1 ) {
        int i_list;
        for ( i_list = 0; i_list <= 4; i_list++ ) {
            double s;          // mandelstam variable "s";
            double t_minus;    // minimum and maximum of the mandelstam variable "t"
            double t_plus;    // as function of the mandelstam variable "s=q2";
            double t_for_s;
            int ijk;        // counter for variable "t";
            int max_ijk;    // maximal value of this counter;

            s = list_of_max_q2_points[i_list];

            t_plus = pow( M1, 2.0 ) + 2.0 * pow( ml, 2.0 ) - s;
            t_plus = t_plus + sqrt( 1.0 - 4.0 * pow( ml, 2.0 ) / s ) *
                                  ( pow( M1, 2.0 ) - s );
            t_plus *= 0.5;

            t_minus = pow( M1, 2.0 ) + 2.0 * pow( ml, 2.0 ) - s;
            t_minus = t_minus - sqrt( 1.0 - 4.0 * pow( ml, 2.0 ) / s ) *
                                    ( pow( M1, 2.0 ) - s );
            t_minus *= 0.5;

            if ( fabs( t_plus - t_minus ) < 0.000001 )
                t_minus = t_plus;

            max_ijk = 1000;
            double dt = ( t_plus - t_minus ) / ( (double)max_ijk );
            if ( fabs( dt ) < 0.00001 )
                dt = 0.0;

            if ( dt < 0.0 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "\n\n In the function EvtbsTollGammaISRFSRAmp::CalcScalarMaxProb(...)"
                    << "\n dt      = " << dt << " < 0."
                    << "\n s       = " << s << "\n t_plus  = " << t_plus
                    << "\n t_minus = " << t_minus << "\n M1      = " << M1
                    << "\n ml      = " << ml << std::endl;
                ::abort();
            }

            for ( ijk = 0; ijk <= max_ijk; ijk++ ) {
                t_for_s = t_minus + dt * ( (double)ijk );

                // B-meson rest frame particles and they kinematics inicialization
                double Eg, El2;
                Eg = ( pow( M1, 2.0 ) - s ) / ( 2.0 * M1 );    // photon energy
                El2 = ( s + t_for_s - pow( ml, 2.0 ) ) /
                      ( 2.0 * M1 );    // ell^- energy

                double modl2;
                modl2 = sqrt( pow( El2, 2.0 ) - pow( ml, 2.0 ) );

                double cosBellminus;    // angle between the B-meson and ell^- directions
                cosBellminus = ( pow( ml, 2.0 ) + 2.0 * Eg * El2 - t_for_s ) /
                               ( 2.0 * Eg * modl2 );

                if ( ( fabs( cosBellminus ) > 1.0 ) &&
                     ( fabs( cosBellminus ) <= 1.0001 ) ) {
                    EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
                        << "\n Debug in the function EvtbsTollGammaISRFSRAmp::CalcMaxProb(...):"
                        << "\n cos(theta) = " << cosBellminus << std::endl;
                    cosBellminus = cosBellminus / fabs( cosBellminus );
                }
                if ( fabs( cosBellminus ) > 1.0001 ) {
                    EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                        << "\n\n In the function EvtbsTollGammaISRFSRAmp::CalcMaxProb(...)"
                        << "\n |cos(theta)| = " << fabs( cosBellminus ) << " > 1"
                        << "\n s       = " << s << "\n t_for_s = " << t_for_s
                        << "\n t_plus  = " << t_plus << "\n t_minus = " << t_minus
                        << "\n dt      = " << dt << "\n Eg      = " << Eg
                        << "\n El2     = " << El2 << "\n modl2   = " << modl2
                        << "\n ml      = " << ml << std::endl;
                    ::abort();
                }

                EvtVector4R p, k, p1, p2;
                p.set( M1, 0.0, 0.0, 0.0 );
                k.set( Eg, Eg, 0.0, 0.0 );
                p2.set( El2, modl2 * cosBellminus,
                        -modl2 * sqrt( 1.0 - pow( cosBellminus, 2.0 ) ), 0.0 );
                p1 = p - k - p2;

                // B-meson state preparation at the rest frame of B-meson
                EvtScalarParticle* scalar_part;
                EvtParticle* root_part;
                scalar_part = new EvtScalarParticle;

                scalar_part->noLifeTime();
                scalar_part->init( parnum, p );
                root_part = (EvtParticle*)scalar_part;
                root_part->setDiagonalSpinDensity();

                // Amplitude initialization
                EvtId listdaug[3];
                listdaug[0] = photnum;
                listdaug[1] = l1num;
                listdaug[2] = l2num;

                EvtAmp amp;
                amp.init( parnum, 3, listdaug );

                // Daughters states preparation at the rest frame of B-meson
                root_part->makeDaughters( 3, listdaug );

                EvtParticle *gamm, *lep1, *lep2;
                gamm = root_part->getDaug( 0 );
                lep1 = root_part->getDaug( 1 );
                lep2 = root_part->getDaug( 2 );

                gamm->noLifeTime();
                lep1->noLifeTime();
                lep2->noLifeTime();

                gamm->init( photnum, k );
                lep1->init( l1num, p1 );
                lep2->init( l2num, p2 );

                EvtSpinDensity rho;
                rho.setDiag( root_part->getSpinStates() );

                // The amplitude calculation at the
                // "maximum amplitude" kinematical configuration
                CalcAmp( root_part, amp, formFactors, WilsCoeff, mu, Nf, sr,
                         res_swch, ias, Egamma_min, CKM_A, CKM_lambda,
                         CKM_barrho, CKM_bareta, mumumass_min );

                // Now find the probability at this q2 and cos theta lepton point
                double nikmax = rho.normalizedProb( amp.getSpinDensity() );

                if ( nikmax > maxfoundprob ) {
                    double maxfoundprob_old;
                    maxfoundprob_old = maxfoundprob;
                    maxfoundprob = nikmax;
                    EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
                        << "\n maxfoundprob ( s =" << s << ",  t = " << t_for_s
                        << " ) = " << maxfoundprob
                        << "\n maxfoundprob_old = " << maxfoundprob_old
                        << "\n ijk =" << ijk << std::endl;
                }

                delete scalar_part;
                //          delete root_part;
                delete gamm;
                delete lep1;
                delete lep2;

            }    // for(ijk=0; ijk<=max_ijk; ijk++)
        }        // i_list - variable loop
    }            // if(res_swch==0||res_swch==1)
    else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n In the function Evtbs2llGammaISRFSRAmp::CalcMaxProb(...)"
            << "\n Unexpected value of the variable res_swch !!!"
            << "\n res_swch = " << res_swch << std::endl;
        ::abort();
    }

    if ( maxfoundprob < 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n In the function Evtbs2llGammaISRFSRAmp::CalcMaxProb(...)"
            << "\n maxfoundprob = " << maxfoundprob << " <0 or =0!"
            << "\n mu =" << mu << " Nf =" << Nf << "\n sr =" << sr
            << " res_swch =" << res_swch << " ias =" << ias
            << "\n Egamma_min =" << Egamma_min << "\n CKM_A      = " << CKM_A
            << " CKM_lambda = " << CKM_lambda << "\n CKM_barrho = " << CKM_barrho
            << " CKM_bareta = " << CKM_bareta << std::endl;
        ::abort();
    }

    if ( maxfoundprob == 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "\n\n In the function Evtbs2llGammaISRFSRAmp::CalcMaxProb(...)"
            << "\n maxfoundprob = " << maxfoundprob << " <0 or =0!"
            << "\n mu =" << mu << " Nf =" << Nf << "\n sr =" << sr
            << " res_swch =" << res_swch << " ias =" << ias
            << "\n Egamma_min =" << Egamma_min << "\n CKM_A      = " << CKM_A
            << " CKM_lambda = " << CKM_lambda << "\n CKM_barrho = " << CKM_barrho
            << " CKM_bareta = " << CKM_bareta << std::endl;
        maxfoundprob = 0.00000001;
    }

    maxfoundprob *= 1.01;

    EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
        << "\n **********************************************************************"
        << "\n The function Evtbs2llGammaISRFSRAmp::CalcMaxProb(...) passed with arguments:"
        << "\n mu =" << mu << " Nf =" << Nf << "\n sr =" << sr
        << " res_swch =" << res_swch << " ias =" << ias
        << "\n CKM_A      = " << CKM_A << " CKM_lambda = " << CKM_lambda
        << "\n Egamma_min =" << Egamma_min << "\n CKM_barrho = " << CKM_barrho
        << " CKM_bareta = " << CKM_bareta
        << "\n The distribution maximum maxfoundprob =" << maxfoundprob
        << "\n **********************************************************************"
        << std::endl;

    return maxfoundprob;
}

// Triangular function
double Evtbs2llGammaISRFSRAmp::lambda( double a, double b, double c )
{
    double l;

    l = pow( a, 2.0 ) + pow( b, 2.0 ) + pow( c, 2.0 ) - 2.0 * a * b -
        2.0 * a * c - 2.0 * b * c;

    return l;
}
