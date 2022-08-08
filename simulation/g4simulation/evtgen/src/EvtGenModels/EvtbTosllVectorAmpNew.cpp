
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

#include "EvtGenModels/EvtbTosllVectorAmpNew.hh"

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

#include "EvtGenModels/EvtbTosllAmpNew.hh"
#include "EvtGenModels/EvtbTosllFFNew.hh"
#include "EvtGenModels/EvtbTosllWilsCoeffNLO.hh"

#include <cstdlib>

//
// The main functiom for the amplitude calculation
//
// input:   *parent      - the pointer to the parent particle (B-meson, the
//                                          object of the EvtParticle class);
//          *formFactors - the pointer to instance of EvtbTosllFFNew class object;
//          *WilsCoeff   - the pointer to the Standart Model Wilson Coefficients class;
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
// return:  amp          - amplitude for the decay  B -> V ell^+ ell^-
//
// Note: in our calculations we assume, that V-meson is the first
//       daughter particle (iV=0) and leptons are the second and thirds
//       daughter particles (il1=1 and il2=2).
//
void EvtbTosllVectorAmpNew::CalcAmp( EvtParticle* parent, EvtAmp& amp,
                                     EvtbTosllFFNew* formFactors,
                                     EvtbTosllWilsCoeffNLO* WilsCoeff,
                                     double mu, int Nf, int res_swch, int ias,
                                     double CKM_A, double CKM_lambda,
                                     double CKM_barrho, double CKM_bareta )
{
    //  FILE *mytest;

    EvtComplex unit1( 1.0, 0.0 );    // real unit
    EvtComplex uniti( 0.0, 1.0 );    // imaginary unit

    int iV = 0;    // V-meson is the first daughter particle
    int il1 = 1,
        il2 = 2;    // leptons are the second and thirds daughter particles

    // transition momentum of the leptonic pair q=k1+k2 or q=p1-p2
    EvtVector4R q = parent->getDaug( il1 )->getP4() +
                    parent->getDaug( il2 )->getP4();

    // Mandelstam variable t=q^2
    double q2 = q.mass2();

    double M1 = parent->mass();                    // B - meson mass
    double M2 = parent->getDaug( iV )->mass();     // V - meson mass
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
    EvtId iddaught = parent->getDaug( iV )->getId();    // The vector meson Id

    // set of the light quark mass value
    if ( ( idparent == EvtPDL::getId( std::string( "B+" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K*+" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B-" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K*-" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K*0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "anti-K*0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "phi" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "phi" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B+" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K_1+" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B-" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K_1-" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K_10" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "anti-K_10" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B+" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K'_1+" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B-" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K'_1-" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K'_10" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "anti-K'_10" ) ) ) ) {
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
           iddaught == EvtPDL::getId( std::string( "rho+" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B-" ) ) &&
           iddaught == EvtPDL::getId( std::string( "rho-" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "rho0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "rho0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "omega" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "omega" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "anti-K*0" ) ) ) ||
         ( idparent == EvtPDL::getId( std::string( "anti-B_s0" ) ) &&
           iddaught == EvtPDL::getId( std::string( "K*0" ) ) ) ) {
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
            << "\n\n The function EvtbTosllVectorAmpNew::CalcAmp(...)"
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

    double a1, a2, a0, v, t1, t2, t3;    // B -> V transition form-factors

    // To get the B -> V transition form-factors
    formFactors->getVectorFF( parent->getId(), parent->getDaug( iV )->getId(),
                              q2, a1, a2, a0, v, t1, t2, t3 );

    // The Wilson Coefficients preparation according to the paper
    // A.J.Buras, M.Munz, Phys.Rev.D52, p.189 (1995)
    EvtComplex c7gam = WilsCoeff->GetC7Eff( mu, Mw, mt, Nf, ias );
    EvtComplex c9eff_b2q = WilsCoeff->GetC9Eff( 0, res_swch, ias, Nf, q2, mb,
                                                ms, mc, mu, mt, Mw, ml,
                                                Relambda_qu, Imlambda_qu );
    EvtComplex c9eff_barb2barq = WilsCoeff->GetC9Eff( 1, res_swch, ias, Nf, q2,
                                                      mb, ms, mc, mu, mt, Mw, ml,
                                                      Relambda_qu, Imlambda_qu );
    EvtComplex c10a = WilsCoeff->GetC10Eff( mt, Mw );

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << "\n\n The function EvtbTosllVectorAmpNew::CalcAmp(...) passed."
    //      << "\n Particle masses:"
    //      << "\n B - meson mass M1 = " << M1
    //      << "\n V - meson mass M2 = " << M2
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
    //      << "   for B -> V transition:"
    //      << "\n v  = " << v
    //      << "\n a0 = " << a0
    //      << "\n a1 = " << a1
    //      << "\n a2 = " << a2
    //      << "\n t1 = " << t1
    //      << "\n t2 = " << t2
    //      << "\n t3 = " << t3
    //      << "\n ============================================================================"
    //      << "\n Wilson Coefficients:"
    //      << "\n Re(c7gam) = " << real(c7gam) << " Im(c7gam) = " << imag(c7gam)
    //      << "\n Re(c9eff_b2q) = " << real(c9eff_b2q)
    //        << " Im(c9eff_b2q) = " << imag(c9eff_b2q)
    //      << "\n Re(c9eff_barb2barq) = " << real(c9eff_barb2barq)
    //        << " Im(c9eff_barb2barq) = " << imag(c9eff_barb2barq)
    //      << "\n Re(c10a)  = " << real(c10a)  << " Im(c10a)  = " << imag(c10a)
    //      << std::endl;

    //  mytest = fopen("output.txt","a");
    //
    //  if(mytest != NULL){
    //     fprintf(mytest,"%lf\n",q2);
    //     fclose(mytest);
    //  }
    //  else{
    //     EvtGenReport(EVTGEN_ERROR,"EvtGen") << "\n Error in writing to file.\n"
    //     << std::endl;
    //     return;
    //  }

    // 4- momentum of the B-meson in the the B-meson rest frame
    EvtVector4R p1 = parent->getP4Restframe();
    EvtVector4R hatp1 = p1 / M1;
    // 4-momentum of the V-meson in the B-meson rest frame
    EvtVector4R p2 = parent->getDaug( 0 )->getP4();
    EvtVector4R hatp2 = p2 / M1;

    double hats = q2 / pow( M1, 2 );
    double hatM2 = M2 / M1;
    double hatmb = mb / M1;
    double hatms = ms / M1;

    // Hadronic matrix element coefficients according to the paper
    // A. Ali, A. Salim Safir, Eur.Phys.J.C25, pp.583-601 (2002)
    // with m_s.NE.0
    EvtComplex a_b2q, a_barb2barq, b_b2q, b_barb2barq, c_b2q, c_barb2barq, e, f,
        g, h;

    a_b2q = 2.0 * c9eff_b2q * v / ( 1.0 + hatM2 ) +
            4.0 * ( hatmb + hatms ) * c7gam * t1 / hats;
    a_barb2barq = 2.0 * c9eff_barb2barq * v / ( 1.0 + hatM2 ) +
                  4.0 * ( hatmb + hatms ) * c7gam * t1 / hats;

    b_b2q = ( c9eff_b2q * a1 +
              2.0 * ( hatmb - hatms ) * ( 1.0 - hatM2 ) * c7gam * t2 / hats ) *
            ( 1.0 + hatM2 );
    b_barb2barq = ( c9eff_barb2barq * a1 + 2.0 * ( hatmb - hatms ) *
                                               ( 1.0 - hatM2 ) * c7gam * t2 /
                                               hats ) *
                  ( 1.0 + hatM2 );

    c_b2q = ( c9eff_b2q * ( 1.0 - hatM2 ) * a2 +
              2.0 * ( hatmb - hatms ) * ( 1.0 - pow( hatM2, 2 ) ) * c7gam * t2 /
                  hats +
              2.0 * ( hatmb - hatms ) * c7gam * t3 ) /
            ( 1 - pow( hatM2, 2 ) );

    c_barb2barq = ( c9eff_barb2barq * ( 1.0 - hatM2 ) * a2 +
                    2.0 * ( hatmb - hatms ) * ( 1.0 - pow( hatM2, 2 ) ) *
                        c7gam * t2 / hats +
                    2.0 * ( hatmb - hatms ) * c7gam * t3 ) /
                  ( 1 - pow( hatM2, 2 ) );

    e = 2.0 * c10a * v / ( 1 + hatM2 );

    f = ( 1.0 + hatM2 ) * c10a * a1;

    g = c10a * a2 / ( 1 + hatM2 );

    h = ( ( 1.0 + hatM2 ) * a1 - ( 1.0 - hatM2 ) * a2 - 2.0 * hatM2 * a0 ) *
        c10a / hats;

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << " a_b2q       = " << a_b2q
    //                          << " a_barb2barq = " << a_barb2barq
    //                          << " b_b2q       = " << b_b2q
    //                          << " b_barb2barq = " << b_barb2barq
    //                          << " c_b2q       = " << c_b2q
    //                          << " c_barb2barq = " << c_barb2barq
    //                          << " e = " << e
    //                          << " f = " << f
    //                          << " g = " << g
    //                          << " h = " << h
    //                          << std::endl;

    // to find ell^+ and ell^- in the B-meson  daughters
    int charge1 = EvtPDL::chg3( parent->getDaug( 1 )->getId() );
    int charge2 = EvtPDL::chg3( parent->getDaug( 2 )->getId() );

    EvtParticle* lepPlus = 0;
    EvtParticle* lepMinus = 0;

    lepPlus = ( charge1 > charge2 ) ? parent->getDaug( 1 ) : parent->getDaug( 2 );
    lepMinus = ( charge1 < charge2 ) ? parent->getDaug( 1 )
                                     : parent->getDaug( 2 );

    EvtTensor4C T1, T2;    // hadronic matrix element tensor structures

    EvtVector4C epsV;    // vector meson polarisation vector

    int i;    // vector meson polarisations counter

    EvtVector4C lvc11, lvc12;    // spin structures for
    EvtVector4C lvc21, lvc22;    // the leptonic vector current

    EvtVector4C lac11, lac12;    // spin structures for
    EvtVector4C lac21, lac22;    // the leptonic axial current

    // B - and barB - mesons descriptors
    EvtIdSet bmesons( "B-", "anti-B0", "anti-B_s0", "B_c-" );
    EvtIdSet bbarmesons( "B+", "B0", "B_s0", "B_c+" );

    EvtId parentID = parent->getId();

    if ( bmesons.contains( parentID ) ) {
        // The amplitude for the decay barB -> barV ell^+ ell^-
        // (b -> q ell^+ ell^- transition)

        T1 = -a_b2q * unit1 * dual( EvtGenFunctions::directProd( hatp1, hatp2 ) ) -
             b_b2q * uniti * EvtTensor4C::g() +
             c_b2q * uniti *
                 EvtGenFunctions::directProd( ( hatp1 + hatp2 ), hatp1 );

        T2 = -e * unit1 * dual( EvtGenFunctions::directProd( hatp1, hatp2 ) ) -
             f * uniti * EvtTensor4C::g() +
             g * uniti * EvtGenFunctions::directProd( ( hatp1 + hatp2 ), hatp1 ) +
             h * uniti * EvtGenFunctions::directProd( ( hatp1 - hatp2 ), hatp1 );

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

        // summing up vector meson polarisations \epsilon^*_{\nu}(i)
        for ( i = 0; i < 3; i++ ) {
            EvtVector4C epsV = parent->getDaug( 0 )->epsParent( i ).conj();

            EvtVector4C E1 = M1 * T1.cont2( epsV );
            EvtVector4C E2 = M1 * T2.cont2( epsV );

            amp.vertex( i, 0, 0, CKM_factor * ( lvc11 * E1 + lac11 * E2 ) );
            amp.vertex( i, 0, 1, CKM_factor * ( lvc12 * E1 + lac12 * E2 ) );
            amp.vertex( i, 1, 0, CKM_factor * ( lvc21 * E1 + lac21 * E2 ) );
            amp.vertex( i, 1, 1, CKM_factor * ( lvc22 * E1 + lac22 * E2 ) );
        }

    } else {
        if ( bbarmesons.contains( parentID ) ) {
            // The amplitude for the decay B -> V ell^+ ell^-
            // (barb -> barq ell^+ ell^- transition)

            T1 = a_barb2barq * unit1 *
                     dual( EvtGenFunctions::directProd( hatp1, hatp2 ) ) -
                 b_barb2barq * uniti * EvtTensor4C::g() +
                 c_barb2barq * uniti *
                     EvtGenFunctions::directProd( ( hatp1 + hatp2 ), hatp1 );

            T2 = e * unit1 * dual( EvtGenFunctions::directProd( hatp1, hatp2 ) ) -
                 f * uniti * EvtTensor4C::g() +
                 g * uniti *
                     EvtGenFunctions::directProd( ( hatp1 + hatp2 ), hatp1 ) +
                 h * uniti *
                     EvtGenFunctions::directProd( ( hatp1 - hatp2 ), hatp1 );

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

            // summing up vector meson polarisations \epsilon^*_{\nu}(i)
            for ( i = 0; i < 3; i++ ) {
                EvtVector4C barepsV = parent->getDaug( 0 )->epsParent( i ).conj();

                EvtVector4C E3 = M1 * T1.cont2( barepsV );
                EvtVector4C E4 = M1 * T2.cont2( barepsV );

                amp.vertex( i, 0, 0,
                            conj( CKM_factor ) * ( lvc11 * E3 + lac11 * E4 ) );
                amp.vertex( i, 0, 1,
                            conj( CKM_factor ) * ( lvc12 * E3 + lac12 * E4 ) );
                amp.vertex( i, 1, 0,
                            conj( CKM_factor ) * ( lvc21 * E3 + lac21 * E4 ) );
                amp.vertex( i, 1, 1,
                            conj( CKM_factor ) * ( lvc22 * E3 + lac22 * E4 ) );
            }

        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "\n\n The function EvtbTosllVectorAmpNew::CalcAmp(...)"
                << "\n Wrong B-meson number" << std::endl;
            ::abort();
        }
    }

    // Test of the signature for Levi-Civita tensor
    //  EvtVector4C Vec0, Vec1, Vec2, Vec3;
    //  EvtTensor4C Ttest;
    //  Vec0.set(1.0,0.0,0.0,0.0);
    //  Vec1.set(0.0,1.0,0.0,0.0);
    //  Vec2.set(0.0,0.0,1.0,0.0);
    //  Vec3.set(0.0,0.0,0.0,1.0);
    //  Ttest=dual(directProd(Vec2,Vec3));
    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << "\n\n\n e^{0123} =" << Ttest.get(0,1) << std::endl;
    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << " e^{1023} =" << Ttest.get(1,0) << std::endl;
    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << " e^{1123} =" << Ttest.get(1,1) << "\n" << std::endl;
    //  EvtVector4C Vtest=Ttest.cont2(Vec1);
    //  for(i=0;i<=3;i++){
    //    EvtGenReport(EVTGEN_NOTICE,"EvtGen") << " Vtest =" << Vtest.get(i) << std::endl;
    //  }
    //  EvtComplex Atest;
    //  Atest=Vec0*Vtest;
    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << "\n Atest =" << Atest << "\n\n\n" << std::endl;
}

//
// The decays B -> V ell^+ ell^- maximum probability calculation for the
// d^2\Gamma/dq^2 d\cos\theta distribution.
//
// \theta - the angle between the vector meson and ell^- directions in the
//          B-meson rest frame.
//
// If ias=0 (nonresonant case), the maximum is achieved at q2=q2_min=4*ml^2.
// If ias=1 (resonant case), the maximum is achieved at q2=M^2_{J/\psi}.
//
double EvtbTosllVectorAmpNew::CalcMaxProb(
    EvtId parnum, EvtId mesnum, EvtId l1num, EvtId l2num,
    EvtbTosllFFNew* formFactors, EvtbTosllWilsCoeffNLO* WilsCoeff, double mu,
    int Nf, int res_swch, int ias, double CKM_A, double CKM_lambda,
    double CKM_barrho, double CKM_bareta )
{
    double maxfoundprob = -100.0;    // maximum of the probability
    int katmax = 0;

    double M1 = EvtPDL::getMeanMass( parnum );    // B - meson mass
    double M2 = EvtPDL::getMeanMass( mesnum );    // V - meson mass
    double ml = EvtPDL::getMeanMass( l1num );     // leptonic mass

    if ( res_swch == 0 ) {
        // B-meson rest frame particles and they kinematics inicialization
        double s_min, t_for_s;
        s_min = 4.0 * pow( ml, 2.0 );
        t_for_s = 0.5 *
                  ( pow( M1, 2.0 ) + pow( M2, 2.0 ) - 2.0 * pow( ml, 2.0 ) );

        double EV, El2;
        EV = ( pow( M1, 2.0 ) + pow( M2, 2.0 ) - s_min ) /
             ( 2.0 * M1 );    // V-meson energy
        El2 = ( s_min + t_for_s - pow( M2, 2.0 ) - pow( ml, 2.0 ) ) /
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
            //       EvtGenReport(EVTGEN_DEBUG,"EvtGen")
            //         << "\n Debug in the function EvtbTosllVectorAmpNew::CalcMaxProb(...):"
            //         << "\n cos(theta) = " << cosVellminus
            //         << std::endl;
            cosVellminus = cosVellminus / fabs( cosVellminus );
        }
        if ( ( modV <= 0.000001 ) || ( modl2 <= 0.000001 ) ) {
            cosVellminus = cosVellminus / fabs( cosVellminus );
            EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
                << "\n Debug in the function EvtbTosllVectorAmpNew::CalcMaxProb(...):"
                << "\n modV       = " << modV << "\n modl2      = " << modl2
                << "\n cos(theta) = " << cosVellminus
                << "\n t_for_s    = " << t_for_s << "\n s_min      = " << s_min
                << "\n EV         = " << EV << "\n El2        = " << El2
                << "\n M2         = " << M2 << "\n ml         = " << ml
                << std::endl;
        }
        if ( fabs( cosVellminus ) > 1.0001 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "\n\n In the function EvtbTosllVectorAmpNew::CalcMaxProb(...)"
                << "\n |cos(theta)| = " << fabs( cosVellminus ) << " > 1"
                << "\n s_min   = " << s_min << "\n t_for_s = " << t_for_s
                << "\n EV      = " << EV << "\n El2     = " << El2
                << "\n modV    = " << modV << "\n modl2   = " << modl2
                << "\n M2      = " << M2 << "\n ml      = " << ml << std::endl;
            ::abort();
        }

        EvtVector4R p1, p2, k1, k2;
        p1.set( M1, 0.0, 0.0, 0.0 );
        p2.set( EV, modV, 0.0, 0.0 );
        k2.set( El2, modl2 * cosVellminus,
                -modl2 * sqrt( 1.0 - pow( cosVellminus, 2.0 ) ), 0.0 );
        k1 = p1 - p2 - k2;

        //     EvtGenReport(EVTGEN_DEBUG,"EvtGen")
        //         << "\n Debug in the function EvtbTosllVectorAmpNew::CalcMaxProb(...):"
        //         << "\n mu =" << mu << " Nf =" << Nf
        //         << " res_swch =" << res_swch
        //         << " ias =" << ias
        //         << "\n CKM_A      = " << CKM_A
        //         << " CKM_lambda = " << CKM_lambda
        //         << "\n CKM_barrho = " << CKM_barrho
        //         << " CKM_bareta = " << CKM_bareta
        //         << "\n M1 = " << M1
        //         << "\n M2 = " << M2
        //         << "\n ml = " << ml
        //         << "\n s_min   = " << s_min
        //         << "\n t_for_s = " << t_for_s
        //         << "\n EV      = " << EV
        //         << "\n El1     = " << El1
        //         << "\n El2     = " << El2
        //         << "\n modV    = " << modV
        //         << "\n modl1   = " << modl1
        //         << "\n modl2   = " << modl2
        //         << "\n cos(theta) = " << cosVellminus
        //         << "\n p1 =" << p1
        //         << "\n p2 =" << p2
        //         << "\n k1 =" << k1
        //         << "\n k2 =" << k2
        //         << std::endl;

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
        CalcAmp( root_part, amp, formFactors, WilsCoeff, mu, Nf, res_swch, ias,
                 CKM_A, CKM_lambda, CKM_barrho, CKM_bareta );

        // Now find the probability at this q2 and cos theta lepton point
        maxfoundprob = rho.normalizedProb( amp.getSpinDensity() );

        delete scalar_part;
        //     delete root_part;
        delete vect;
        delete lep1;
        delete lep2;

    }    // if(res_swch==0)

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
        for ( k = 0; k <= 1000; k++ ) {
            t_for_s = t_plus - dt * ( (double)k );

            if ( ( t_for_s < t_minus ) && ( t_for_s >= ( 0.9999 * t_minus ) ) ) {
                t_for_s = t_minus;
            }
            if ( t_for_s < ( 0.9999 * t_minus ) ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "\n\n In the function EvtbTosllVectorAmpNew::CalcMaxProb(...)"
                    << "\n t_for_s = " << t_for_s << " < t_minus = " << t_minus
                    << " ! "
                    << "\n t_plus  = " << t_plus << "\n dt      = " << dt
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
                //         EvtGenReport(EVTGEN_DEBUG,"EvtGen")
                //           << "\n Debug in the function EvtbTosllVectorAmpNew::CalcMaxProb(...):"
                //           << "\n cos(theta) = " << cosVellminus
                //           << std::endl;
                cosVellminus = cosVellminus / fabs( cosVellminus );
            }
            if ( ( modV <= 0.000001 ) || ( modl2 <= 0.000001 ) ) {
                cosVellminus = cosVellminus / fabs( cosVellminus );
                EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
                    << "\n Debug in the function EvtbTosllVectorAmpNew::CalcMaxProb(...):"
                    << "\n modV       = " << modV << "\n modl2      = " << modl2
                    << "\n cos(theta) = " << cosVellminus
                    << "\n s          = " << s << "\n t_for_s    = " << t_for_s
                    << "\n t_plus     = " << t_plus
                    << "\n t_minus    = " << t_minus << "\n dt         = " << dt
                    << "\n EV         = " << EV << "\n El2        = " << El2
                    << "\n M2         = " << M2 << "\n ml         = " << ml
                    << std::endl;
            }
            if ( fabs( cosVellminus ) > 1.0001 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "\n\n In the function EvtbTosllVectorAmpNew::CalcMaxProb(...)"
                    << "\n |cos(theta)| = " << fabs( cosVellminus ) << " > 1"
                    << "\n s       = " << s << "\n t_for_s = " << t_for_s
                    << "\n EV      = " << EV << "\n El2     = " << El2
                    << "\n modV    = " << modV << "\n modl2   = " << modl2
                    << "\n M2      = " << M2 << "\n ml      = " << ml
                    << std::endl;
                ::abort();
            }

            EvtVector4R p1, p2, k1, k2;
            p1.set( M1, 0.0, 0.0, 0.0 );
            p2.set( EV, modV, 0.0, 0.0 );
            k2.set( El2, modl2 * cosVellminus,
                    -modl2 * sqrt( 1.0 - pow( cosVellminus, 2.0 ) ), 0.0 );
            k1 = p1 - p2 - k2;

            //       EvtGenReport(EVTGEN_DEBUG,"EvtGen")
            //           << "\n Debug in the function EvtbTosllVectorAmpNew::CalcMaxProb(...):"
            //           << "\n mu =" << mu << " Nf =" << Nf
            //           << " res_swch =" << res_swch
            //           << " ias =" << ias
            //           << "\n CKM_A      = " << CKM_A
            //           << " CKM_lambda = " << CKM_lambda
            //           << "\n CKM_barrho = " << CKM_barrho
            //           << " CKM_bareta = " << CKM_bareta
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
                     ias, CKM_A, CKM_lambda, CKM_barrho, CKM_bareta );

            // Now find the probability at this q2 and cos theta lepton point
            double nikmax = rho.normalizedProb( amp.getSpinDensity() );

            if ( nikmax > maxfoundprob ) {
                maxfoundprob = nikmax;
                katmax = k;
                EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
                    << "\n maxfoundprob ( s =" << s << ",  t = " << t_for_s
                    << " ) = " << maxfoundprob << "\n k =" << katmax
                    << std::endl;
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
            << "\n\n In the function EvtbTosllVectorAmpNew::CalcMaxProb(...)"
            << "\n maxfoundprob = " << maxfoundprob << " <0 or =0!"
            << "\n res_swch     = " << res_swch << std::endl;
        ::abort();
    }

    EvtGenReport( EVTGEN_NOTICE, "EvtGen" )
        << "\n maxfoundprob (...) = " << maxfoundprob << std::endl;

    maxfoundprob *= 1.01;

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen")
    //         << "\n ***************************************************************************"
    //         << "\n The function EvtbTosllVectorAmpNew::CalcMaxProb(...) passed with arguments:"
    //         << "\n mu =" << mu << " Nf =" << Nf
    //         << " res_swch =" << res_swch
    //         << " ias =" << ias
    //         << "\n CKM_A      = " << CKM_A
    //         << " CKM_lambda = " << CKM_lambda
    //         << "\n CKM_barrho = " << CKM_barrho
    //         << " CKM_bareta = " << CKM_bareta
    //         << "\n The distribution maximum maxfoundprob =" << maxfoundprob
    //         << "\n k = " << katmax
    //         << "\n ***************************************************************************"
    //         << std::endl;

    return maxfoundprob;
}

// Triangular function
double EvtbTosllVectorAmpNew::lambda( double a, double b, double c )
{
    double l;

    l = pow( a, 2.0 ) + pow( b, 2.0 ) + pow( c, 2.0 ) - 2.0 * a * b -
        2.0 * a * c - 2.0 * b * c;

    return l;
}
