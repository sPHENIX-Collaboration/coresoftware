
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

#include "EvtGenModels/EvtBsMuMuKK.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtdFunction.hh"

const double pi = EvtConst::pi;
const EvtComplex I = EvtComplex( 0.0, 1.0 );
const double sq2 = sqrt( 2.0 );

std::string EvtBsMuMuKK::getName()
{
    return "BS_MUMUKK";
}

EvtDecayBase* EvtBsMuMuKK::clone()
{
    return new EvtBsMuMuKK;
}

void EvtBsMuMuKK::init()
{
    // DecFile parameters
    checkNArg( 37 );

    // Non-resonant S wave
    f_S_NR = getArg( 0 );
    delta_S_NR = getArg( 1 );
    phis_S_NR = getArg( 2 );
    lambda_S_NR_abs = getArg( 3 );

    // f0 (S wave)
    f_f0 = getArg( 4 );
    delta_f0 = getArg( 5 );
    phis_f0 = getArg( 6 );
    lambda_f0_abs = getArg( 7 );

    // phi (P wave)
    f_phi = getArg( 8 );
    f_phi_0 = getArg( 9 );
    delta_phi_0 = getArg( 10 );
    phis_phi_0 = getArg( 11 );
    lambda_phi_0_abs = getArg( 12 );
    f_phi_perp = getArg( 13 );
    delta_phi_perp = pi - getArg( 14 );
    phis_phi_perp = getArg( 15 );
    lambda_phi_perp_abs = getArg( 16 );
    delta_phi_par = pi - getArg( 17 );
    phis_phi_par = getArg( 18 );
    lambda_phi_par_abs = getArg( 19 );

    // f2' (D wave)
    f_f2p_0 = getArg( 20 );
    delta_f2p_0 = getArg( 21 );
    phis_f2p_0 = getArg( 22 );
    lambda_f2p_0_abs = getArg( 23 );
    f_f2p_perp = getArg( 24 );
    delta_f2p_perp = pi - getArg( 25 );
    phis_f2p_perp = getArg( 26 );
    lambda_f2p_perp_abs = getArg( 27 );
    delta_f2p_par = pi - getArg( 28 );
    phis_f2p_par = getArg( 29 );
    lambda_f2p_par_abs = getArg( 30 );

    // Time dependence
    Gamma = getArg( 31 );
    deltaGamma = getArg( 32 );
    deltaMs = getArg( 33 );

    // mKK window
    Mf0 = getArg( 34 );
    kin_lower_limit = getArg( 35 );    // the minimum is approx 2.03*MKp
    kin_upper_limit = getArg( 36 );

    // PDG masses
    MBs = EvtPDL::getMass( EvtPDL::getId( "B_s0" ) );
    MJpsi = EvtPDL::getMeanMass( EvtPDL::getId( "J/psi" ) );
    Mphi = EvtPDL::getMeanMass( EvtPDL::getId( "phi" ) );
    Mf2p = EvtPDL::getMeanMass( EvtPDL::getId( "f'_2" ) );
    MKp = EvtPDL::getMass( EvtPDL::getId( "K+" ) );
    MKm = EvtPDL::getMass( EvtPDL::getId( "K-" ) );
    MK0 = EvtPDL::getMass( EvtPDL::getId( "K0" ) );
    Mpip = EvtPDL::getMass( EvtPDL::getId( "pi+" ) );
    Mpi0 = EvtPDL::getMass( EvtPDL::getId( "pi0" ) );
    Mmu = EvtPDL::getMass( EvtPDL::getId( "mu+" ) );

    double MBsSq = MBs * MBs;

    // Amplitudes and other time parameters
    A_S_NR = sqrt( f_S_NR );
    A_f0 = sqrt( f_f0 );

    A_phi_0 = sqrt( f_phi_0 * f_phi );
    A_phi_perp = sqrt( f_phi_perp * f_phi );
    // Use fabs to make sure subtractions are >= 0, since subtracting 0 from 0 can give -0
    A_phi_par = sqrt(
        fabs( f_phi - A_phi_perp * A_phi_perp - A_phi_0 * A_phi_0 ) );

    f_f2p = fabs( 1.0 - f_S_NR - f_f0 - f_phi );
    A_f2p_0 = sqrt( f_f2p_0 * f_f2p );
    A_f2p_perp = sqrt( f_f2p_perp * f_f2p );
    A_f2p_par = sqrt(
        fabs( f_f2p - A_f2p_perp * A_f2p_perp - A_f2p_0 * A_f2p_0 ) );

    ctau = 1.0 / Gamma;
    Gamma0phi = EvtPDL::getWidth( EvtPDL::getId( "phi" ) );
    Gamma0f2p = EvtPDL::getWidth( EvtPDL::getId( "f'_2" ) );

    kin_middle = 0.5 * ( kin_upper_limit + kin_lower_limit );

    int_const_NR = sqrt(
        Integral( 1.0, 1.0, 0, 1, 1.0, kin_lower_limit, kin_upper_limit, 0 ) );

    int_Flatte_f0 = sqrt(
        Integral( 1.0, Mf0, 0, 1, 1.0, kin_lower_limit, kin_upper_limit, 1 ) );

    p30Kp_mid_CMS = sqrt( ( pow( kin_middle, 2 ) - pow( MKp + MKm, 2 ) ) *
                          ( pow( kin_middle, 2 ) - pow( MKp - MKm, 2 ) ) ) /
                    ( 2.0 * kin_middle );

    p30Kp_ll_CMS = sqrt( ( pow( kin_lower_limit, 2 ) - pow( MKp + MKm, 2 ) ) *
                         ( pow( kin_lower_limit, 2 ) - pow( MKp - MKm, 2 ) ) ) /
                   ( 2.0 * kin_lower_limit );

    p30Kp_phi_CMS = sqrt( ( Mphi * Mphi - pow( MKp + MKm, 2 ) ) *
                          ( Mphi * Mphi - pow( MKp - MKm, 2 ) ) ) /
                    ( 2.0 * Mphi );

    p30Kp_f2p_CMS = sqrt( ( Mf2p * Mf2p - pow( MKp + MKm, 2 ) ) *
                          ( Mf2p * Mf2p - pow( MKp - MKm, 2 ) ) ) /
                    ( 2.0 * Mf2p );

    p30Jpsi_mid_CMS = sqrt( ( MBsSq - pow( kin_middle + MJpsi, 2 ) ) *
                            ( MBsSq - pow( kin_middle - MJpsi, 2 ) ) ) /
                      ( 2.0 * MBs );

    p30Jpsi_ll_CMS = sqrt( ( MBsSq - pow( kin_lower_limit + MJpsi, 2 ) ) *
                           ( MBsSq - pow( kin_lower_limit - MJpsi, 2 ) ) ) /
                     ( 2.0 * MBs );

    p30Jpsi_phi_CMS = sqrt( ( MBsSq - pow( Mphi + MJpsi, 2 ) ) *
                            ( MBsSq - pow( Mphi - MJpsi, 2 ) ) ) /
                      ( 2.0 * MBs );

    p30Jpsi_f2p_CMS = sqrt( ( MBsSq - pow( Mf2p + MJpsi, 2 ) ) *
                            ( MBsSq - pow( Mf2p - MJpsi, 2 ) ) ) /
                      ( 2.0 * MBs );

    int_BW_phi = sqrt( Integral( Gamma0phi, Mphi, 1, 0, p30Kp_phi_CMS,
                                 kin_lower_limit, kin_upper_limit, 2 ) );

    int_BW_f2p = sqrt( Integral( Gamma0f2p, Mf2p, 2, 1, p30Kp_f2p_CMS,
                                 kin_lower_limit, kin_upper_limit, 2 ) );

    // 4 daughters
    checkNDaug( 4 );

    // Spin-0 parent
    checkSpinParent( EvtSpinType::SCALAR );    // B_s0 (anti-B_s0)

    // Daughters
    checkSpinDaughter( 0, EvtSpinType::DIRAC );     // mu+  (mu-)
    checkSpinDaughter( 1, EvtSpinType::DIRAC );     // mu-  (mu+)
    checkSpinDaughter( 2, EvtSpinType::SCALAR );    // K+   (K-)
    checkSpinDaughter( 3, EvtSpinType::SCALAR );    // K-   (K+)

    // B_s0 parent (Parent must be B_s0 or anti-B_s0)
    const EvtId p = getParentId();
    if ( p != EvtPDL::getId( "B_s0" ) && p != EvtPDL::getId( "anti-B_s0" ) ) {
        assert( 0 );
    }

    // Daughter types and ordering (should be mu+-, mu-+, K+-, K-+)
    const EvtId d1 = getDaug( 0 );
    const EvtId d2 = getDaug( 1 );
    const EvtId d3 = getDaug( 2 );
    const EvtId d4 = getDaug( 3 );
    if ( !( ( d1 == EvtPDL::getId( "mu+" ) || d1 == EvtPDL::getId( "mu-" ) ) &&
            ( d2 == EvtPDL::getId( "mu-" ) || d2 == EvtPDL::getId( "mu+" ) ) &&
            ( d3 == EvtPDL::getId( "K+" ) || d3 == EvtPDL::getId( "K-" ) ) &&
            ( d4 == EvtPDL::getId( "K-" ) || d4 == EvtPDL::getId( "K+" ) ) ) ) {
        assert( 0 );
    }
}

// Get ProbMax
void EvtBsMuMuKK::initProbMax()
{
    const EvtComplex term11 = sqrt( p30Jpsi_f2p_CMS * p30Kp_f2p_CMS );

    const EvtComplex term12 = X_J( 2, p30Kp_f2p_CMS, 0 ) *
                              X_J( 1, p30Jpsi_f2p_CMS, 1 ) * p30Kp_f2p_CMS *
                              p30Kp_f2p_CMS * p30Jpsi_f2p_CMS *
                              ( A_f2p_0 + 0.3 * A_f2p_perp + 0.3 * A_f2p_par );

    const EvtComplex term13 = f_f2p *
                              Breit_Wigner( Gamma0f2p, Mf2p, Mf2p, 2,
                                            p30Kp_f2p_CMS, p30Kp_f2p_CMS ) /
                              int_BW_f2p;

    const EvtComplex term21 = sqrt( p30Jpsi_phi_CMS * p30Kp_phi_CMS );

    const EvtComplex term22 = X_J( 1, p30Kp_phi_CMS, 0 ) * p30Kp_phi_CMS *
                              ( 0.65 * A_phi_0 + 0.6 * A_phi_perp +
                                0.6 * A_phi_par );

    const EvtComplex term23 = f_phi *
                              Breit_Wigner( Gamma0phi, Mphi, Mphi, 1,
                                            p30Kp_phi_CMS, p30Kp_phi_CMS ) /
                              int_BW_phi;

    const EvtComplex term31 = sqrt( p30Jpsi_ll_CMS * p30Kp_ll_CMS );

    const EvtComplex term32 = X_J( 1, p30Jpsi_ll_CMS, 1 ) * p30Jpsi_ll_CMS;

    const EvtComplex term33 = f_f0 * Flatte( Mf0, kin_lower_limit ) /
                              int_Flatte_f0;

    const EvtComplex term41 = sqrt( p30Jpsi_mid_CMS * p30Kp_mid_CMS );

    const EvtComplex term42 = X_J( 1, p30Jpsi_mid_CMS, 1 ) * p30Jpsi_mid_CMS;

    const EvtComplex term43 = 1.2 * f_S_NR / int_const_NR;

    const EvtComplex hm = term11 * term12 * term13 + term21 * term22 * term23 +
                          term31 * term32 * term33 + term41 * term42 * term43;

    // Increase by 10%
    setProbMax( 0.5 * abs2( hm ) * 1.1 );
}

// Decay function
void EvtBsMuMuKK::decay( EvtParticle* p )
{
    EvtId other_b;
    double time( 0.0 );
    EvtCPUtil::getInstance()->OtherB( p, time, other_b );
    time = -log( EvtRandom::Flat() ) *
           ctau;    // This overrules the ctau made in OtherB

    if ( EvtCPUtil::getInstance()->isBsMixed( p ) ) {
        p->getParent()->setLifetime( time * EvtConst::c / 1e12 );    // units: mm
    } else {
        p->setLifetime( time * EvtConst::c / 1e12 );    // units: mm
    }

    double DGtime = 0.25 * deltaGamma * time;
    double DMtime = 0.5 * deltaMs * time;
    double mt = exp( -DGtime );
    double pt = exp( +DGtime );
    double cDMt = cos( DMtime );
    double sDMt = sin( DMtime );
    EvtComplex termplus = EvtComplex( cDMt, sDMt );
    EvtComplex terminus = EvtComplex( cDMt, -sDMt );

    EvtComplex gplus = 0.5 * ( mt * termplus + pt * terminus );
    EvtComplex gminus = 0.5 * ( mt * termplus - pt * terminus );

    EvtId BSB = EvtPDL::getId( "anti-B_s0" );

    // Flavour: first assume B_s0, otherwise choose anti-B_s0
    int q( 1 );
    if ( other_b == BSB ) {
        q = -1;
    }
    p->setAttribute( "q", q );

    // Amplitudes
    EvtComplex a_S_NR = AmpTime( q, gplus, gminus, delta_S_NR, lambda_S_NR_abs,
                                 A_S_NR, phis_S_NR, -1 );

    EvtComplex a_f0 = AmpTime( q, gplus, gminus, delta_f0, lambda_f0_abs, A_f0,
                               phis_f0, -1 );

    EvtComplex a0_phi = AmpTime( q, gplus, gminus, delta_phi_0,
                                 lambda_phi_0_abs, A_phi_0, phis_phi_0, 1 );

    EvtComplex aperp_phi = AmpTime( q, gplus, gminus, delta_phi_perp,
                                    lambda_phi_perp_abs, A_phi_perp,
                                    phis_phi_perp, -1 );

    EvtComplex apar_phi = AmpTime( q, gplus, gminus, delta_phi_par,
                                   lambda_phi_par_abs, A_phi_par, phis_phi_par,
                                   1 );

    EvtComplex a0_f2p = AmpTime( q, gplus, gminus, delta_f2p_0,
                                 lambda_f2p_0_abs, A_f2p_0, phis_f2p_0, -1 );

    EvtComplex aperp_f2p = AmpTime( q, gplus, gminus, delta_f2p_perp,
                                    lambda_f2p_perp_abs, A_f2p_perp,
                                    phis_f2p_perp, 1 );

    EvtComplex apar_f2p = AmpTime( q, gplus, gminus, delta_f2p_par,
                                   lambda_f2p_par_abs, A_f2p_par, phis_f2p_par,
                                   -1 );

    // Generate 4-momenta
    double mKK = EvtRandom::Flat( kin_lower_limit, kin_upper_limit );
    double mass[10] = {MJpsi, mKK, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double Kmass[10] = {MKp, MKm, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double muMass[10] = {Mmu, Mmu, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    EvtVector4R mypV[2], mypK[2], mypmu[2];
    EvtGenKine::PhaseSpace( 2, mass, mypV, MBs );
    EvtGenKine::PhaseSpace( 2, Kmass, mypK, mKK );
    EvtGenKine::PhaseSpace( 2, muMass, mypmu, MJpsi );

    EvtVector4R p4mup = boostTo( mypmu[0], mypV[0] );
    EvtVector4R p4mum = boostTo( mypmu[1], mypV[0] );
    EvtVector4R p4Kp = boostTo( mypK[0], mypV[1] );
    EvtVector4R p4Km = boostTo( mypK[1], mypV[1] );

    p->makeDaughters( getNDaug(), getDaugs() );

    EvtParticle* thisparticle;
    EvtParticle *muplus, *muminus, *Kplus, *Kminus;

    // Check particle ID
    for ( int k = 0; k <= 3; k++ ) {
        thisparticle = p->getDaug( k );
        EvtId pId = thisparticle->getId();

        if ( pId == EvtPDL::getId( "mu+" ) ) {
            muplus = thisparticle;
            muplus->init( getDaug( k ), p4mup );

        } else if ( pId == EvtPDL::getId( "mu-" ) ) {
            muminus = thisparticle;
            muminus->init( getDaug( k ), p4mum );

        } else if ( pId == EvtPDL::getId( "K+" ) ) {
            Kplus = thisparticle;
            Kplus->init( getDaug( k ), p4Kp );

        } else if ( pId == EvtPDL::getId( "K-" ) ) {
            Kminus = thisparticle;
            Kminus->init( getDaug( k ), p4Km );
        }
    }

    EvtVector4R p4KK = p4Kp + p4Km;
    EvtVector4R p4mumu = p4mup + p4mum;
    EvtVector4R p4Bs = p4mumu + p4KK;

    double p4KK_mass2 = p4KK.mass2();
    double p4KK_mass = p4KK.mass();
    double p4Bs_mass2 = p4Bs.mass2();
    double p4Bs_mass = p4Bs.mass();

    // Kp momentum in the KK CMS
    double p3Kp_KK_CMS = sqrt( ( p4KK_mass2 - pow( MKp + MKm, 2 ) ) *
                               ( p4KK_mass2 - pow( MKp - MKm, 2 ) ) ) /
                         ( 2.0 * p4KK_mass );

    // J/psi momentum in the KK CMS
    double p3Jpsi_KK_CMS = sqrt( ( p4Bs_mass2 - pow( p4KK_mass + MJpsi, 2 ) ) *
                                 ( p4Bs_mass2 - pow( p4KK_mass - MJpsi, 2 ) ) ) /
                           ( 2.0 * p4Bs_mass );

    // Mass lineshapes

    // Non-resonant S wave
    EvtComplex P_NR = 1.0 / int_const_NR;

    // f0 Flatte
    EvtComplex F_f0 = Flatte( Mf0, p4KK_mass ) / int_Flatte_f0;

    // phi Breit Wigner
    EvtComplex BW_phi = Breit_Wigner( Gamma0phi, Mphi, p4KK_mass, 1,
                                      p30Kp_phi_CMS, p3Kp_KK_CMS ) /
                        int_BW_phi;

    // f2' Breit Wigner
    EvtComplex BW_f2p = Breit_Wigner( Gamma0f2p, Mf2p, p4KK_mass, 1,
                                      p30Kp_f2p_CMS, p3Kp_KK_CMS ) /
                        int_BW_f2p;

    // Barrier factors: Always taking the lowest Bs L
    double X_KK_0 = 1.0;
    double X_KK_1 = X_J( 1, p3Kp_KK_CMS, 0 );
    double X_KK_2 = X_J( 2, p3Kp_KK_CMS, 0 );
    double X_NR_Jpsi_1 = X_J( 1, p3Jpsi_KK_CMS, 1 );
    double X_f0_Jpsi_1 = X_J( 1, p3Jpsi_KK_CMS, 1 );
    double X_phi_Jpsi_0 = 1.0;
    double X_f2p_Jpsi_1 = X_J( 1, p3Jpsi_KK_CMS, 1 );

    // Birth momentum factors: pow(p3(K+),LR)* pow(p3(J/psi),LB)
    double f_PHSP = sqrt( p3Jpsi_KK_CMS * p3Kp_KK_CMS );
    double f_BMF_NR = p3Jpsi_KK_CMS;
    double f_BMF_f0 = p3Jpsi_KK_CMS;
    double f_BMF_phi = p3Kp_KK_CMS;
    double f_BMF_f2p = p3Kp_KK_CMS * p3Kp_KK_CMS * p3Jpsi_KK_CMS;

    // Angular distribution and sum over KK states
    double CosK = EvtDecayAngle( p4Bs, p4KK, p4Kp );
    double CosMu = EvtDecayAngle( p4Bs, p4mumu, p4mup );
    double chi = EvtDecayAngleChi( p4Bs, p4mup, p4mum, p4Kp, p4Km );

    // Build helicity amplitudes

    // phi
    EvtComplex H0_phi = a0_phi;
    EvtComplex Hp_phi = ( apar_phi + aperp_phi ) / sq2;
    EvtComplex Hm_phi = ( apar_phi - aperp_phi ) / sq2;

    // f2p
    EvtComplex H0_f2p = a0_f2p;
    EvtComplex Hp_f2p = ( apar_f2p + aperp_f2p ) / sq2;
    EvtComplex Hm_f2p = ( apar_f2p - aperp_f2p ) / sq2;

    // muon polarization +1
    EvtComplex swaveangdist1 = AngularDist( 0, 0, 1, CosK, CosMu, chi );

    // KK Spin-0 NR
    EvtComplex mp_hS_NR = a_S_NR * swaveangdist1;
    EvtComplex Amp_p_NR = P_NR * X_KK_0 * X_NR_Jpsi_1 * f_BMF_NR * mp_hS_NR;

    // KK Spin-0 f0
    EvtComplex mp_h_f0 = a_f0 * swaveangdist1;
    EvtComplex Amp_p_f0 = F_f0 * X_KK_0 * X_f0_Jpsi_1 * f_BMF_f0 * mp_h_f0;

    // KK Spin-1
    EvtComplex mp_h0_phi = H0_phi * AngularDist( 1, 0, 1, CosK, CosMu, chi );
    EvtComplex mp_hp_phi = Hp_phi * AngularDist( 1, 1, 1, CosK, CosMu, chi );
    EvtComplex mp_hm_phi = Hm_phi * AngularDist( 1, -1, 1, CosK, CosMu, chi );
    EvtComplex Amp_p_phi = BW_phi * X_KK_1 * X_phi_Jpsi_0 * f_BMF_phi *
                           ( mp_h0_phi + mp_hp_phi + mp_hm_phi );

    // KK Spin-2
    EvtComplex mp_h0_f2p = H0_f2p * AngularDist( 2, 0, 1, CosK, CosMu, chi );
    EvtComplex mp_hp_f2p = Hp_f2p * AngularDist( 2, 1, 1, CosK, CosMu, chi );
    EvtComplex mp_hm_f2p = Hm_f2p * AngularDist( 2, -1, 1, CosK, CosMu, chi );
    EvtComplex Amp_p_f2p = BW_f2p * X_KK_2 * X_f2p_Jpsi_1 * f_BMF_f2p *
                           ( mp_h0_f2p + mp_hp_f2p + mp_hm_f2p );

    // muon polarization -1
    EvtComplex swaveangdist2 = AngularDist( 0, 0, -1, CosK, CosMu, chi );

    // KK Spin-0 NR
    EvtComplex mm_hS_NR = a_S_NR * swaveangdist2;
    EvtComplex Amp_m_NR = P_NR * X_KK_0 * X_NR_Jpsi_1 * f_BMF_NR * mm_hS_NR;

    // KK Spin-0
    EvtComplex mm_h_f0 = a_f0 * swaveangdist2;
    EvtComplex Amp_m_f0 = F_f0 * X_KK_0 * X_f0_Jpsi_1 * f_BMF_f0 * mm_h_f0;

    // KK Spin-1
    EvtComplex mm_h0_phi = H0_phi * AngularDist( 1, 0, -1, CosK, CosMu, chi );
    EvtComplex mm_hp_phi = Hp_phi * AngularDist( 1, +1, -1, CosK, CosMu, chi );
    EvtComplex mm_hm_phi = Hm_phi * AngularDist( 1, -1, -1, CosK, CosMu, chi );
    EvtComplex Amp_m_phi = BW_phi * X_KK_1 * X_phi_Jpsi_0 * f_BMF_phi *
                           ( mm_h0_phi + mm_hp_phi + mm_hm_phi );

    // KK Spin-2
    EvtComplex mm_h0_f2p = H0_f2p * AngularDist( 2, 0, -1, CosK, CosMu, chi );
    EvtComplex mm_hp_f2p = Hp_f2p * AngularDist( 2, 1, -1, CosK, CosMu, chi );
    EvtComplex mm_hm_f2p = Hm_f2p * AngularDist( 2, -1, -1, CosK, CosMu, chi );
    EvtComplex Amp_m_f2p = BW_f2p * X_KK_2 * X_f2p_Jpsi_1 * f_BMF_f2p *
                           ( mm_h0_f2p + mm_hp_f2p + mm_hm_f2p );

    // Total amplitudes
    EvtComplex Amp_tot_plus = f_PHSP *
                              ( Amp_p_NR + Amp_p_f0 + Amp_p_phi + Amp_p_f2p );
    EvtComplex Amp_tot_minus = f_PHSP *
                               ( Amp_m_NR + Amp_m_f0 + Amp_m_phi + Amp_m_f2p );

    vertex( 0, 0, 0.0 );
    vertex( 0, 1, Amp_tot_plus );
    vertex( 1, 0, Amp_tot_minus );
    vertex( 1, 1, 0.0 );
}

// Rho function
EvtComplex EvtBsMuMuKK::GetRho( const double m0, const double m ) const
{
    double rho_sq = 1.0 - ( 4.0 * m0 * m0 / ( m * m ) );
    EvtComplex rho;

    if ( rho_sq > 0.0 ) {
        rho = EvtComplex( sqrt( rho_sq ), 0.0 );
    } else {
        rho = EvtComplex( 0.0, sqrt( -rho_sq ) );
    }

    return rho;
}

// Flatte function
EvtComplex EvtBsMuMuKK::Flatte( const double m0, const double m ) const
{
    double gpipi = 0.167;
    double gKK = 3.05 * gpipi;

    EvtComplex term1 = ( 2.0 * GetRho( Mpip, m ) + GetRho( Mpi0, m ) ) / 3.0;
    EvtComplex term2 = ( GetRho( MKp, m ) + GetRho( MK0, m ) ) / 2.0;

    EvtComplex w = gpipi * term1 + gKK * term2;

    EvtComplex Flatte_0 = 1.0 / ( m0 * m0 - m * m - I * m0 * w );

    return Flatte_0;
}

// Breit-Wigner function
EvtComplex EvtBsMuMuKK::Breit_Wigner( const double Gamma0, const double m0,
                                      const double m, const int J,
                                      const double q0, const double q ) const
{
    double X_J_q0_sq = pow( X_J( J, q0, 0 ), 2 );
    double X_J_q_sq = pow( X_J( J, q, 0 ), 2 );

    double Gamma = Gamma0 * pow( q / q0, 2 * J + 1 ) * ( m0 / m ) *
                   ( X_J_q_sq / X_J_q0_sq );

    return 1.0 / ( m0 * m0 - m * m - I * m0 * Gamma );
}

// Integral
double EvtBsMuMuKK::Integral( const double Gamma0, const double m0, const int JR,
                              const int JB, const double q0, const double M_KK_ll,
                              const double M_KK_ul, const int fcntype ) const
{
    int bins = 1000;
    double bin_width = ( M_KK_ul - M_KK_ll ) / static_cast<double>( bins );
    EvtComplex integral( 0.0, 0.0 );
    double sumMKpKm2 = pow( MKp + MKm, 2 );
    double diffMKpKm2 = pow( MKp - MKm, 2 );
    double MBs2 = pow( MBs, 2 );

    for ( int i = 0; i < bins; i++ ) {
        double M_KK_i = M_KK_ll + static_cast<double>( i ) * bin_width;
        double M_KK_f = M_KK_ll + static_cast<double>( i + 1 ) * bin_width;
        double M_KK_i_sq = M_KK_i * M_KK_i;
        double M_KK_f_sq = M_KK_f * M_KK_f;

        double p3Kp_KK_CMS_i = sqrt( ( M_KK_i_sq - sumMKpKm2 ) *
                                     ( M_KK_i_sq - diffMKpKm2 ) ) /
                               ( 2.0 * M_KK_i );
        double p3Kp_KK_CMS_f = sqrt( ( M_KK_f_sq - sumMKpKm2 ) *
                                     ( M_KK_f_sq - diffMKpKm2 ) ) /
                               ( 2.0 * M_KK_f );

        double p3Jpsi_Bs_CMS_i = sqrt( ( MBs2 - pow( M_KK_i + MJpsi, 2 ) ) *
                                       ( MBs2 - pow( M_KK_i - MJpsi, 2 ) ) ) /
                                 ( 2.0 * MBs );
        double p3Jpsi_Bs_CMS_f = sqrt( ( MBs2 - pow( M_KK_f + MJpsi, 2 ) ) *
                                       ( MBs2 - pow( M_KK_f - MJpsi, 2 ) ) ) /
                                 ( 2.0 * MBs );

        double f_PHSP_i = sqrt( p3Kp_KK_CMS_i * p3Jpsi_Bs_CMS_i );
        double f_PHSP_f = sqrt( p3Kp_KK_CMS_f * p3Jpsi_Bs_CMS_f );

        double f_MBF_KK_i = pow( p3Kp_KK_CMS_i, JR );
        double f_MBF_KK_f = pow( p3Kp_KK_CMS_f, JR );

        double f_MBF_Bs_i = pow( p3Jpsi_Bs_CMS_i, JB );
        double f_MBF_Bs_f = pow( p3Jpsi_Bs_CMS_f, JB );

        double X_JR_i = X_J( JR, p3Kp_KK_CMS_i, 0 );
        double X_JR_f = X_J( JR, p3Kp_KK_CMS_f, 0 );

        double X_JB_i = X_J( JB, p3Jpsi_Bs_CMS_i, 1 );
        double X_JB_f = X_J( JB, p3Jpsi_Bs_CMS_f, 1 );

        EvtComplex fcn_i( 1.0, 0.0 ), fcn_f( 1.0, 0.0 );

        if ( fcntype == 1 ) {
            fcn_i = Flatte( m0, M_KK_i );
            fcn_f = Flatte( m0, M_KK_f );

        } else if ( fcntype == 2 ) {
            fcn_i = Breit_Wigner( Gamma0, m0, M_KK_i, JR, q0, p3Kp_KK_CMS_i );
            fcn_f = Breit_Wigner( Gamma0, m0, M_KK_f, JR, q0, p3Kp_KK_CMS_f );
        }

        EvtComplex a_i = f_PHSP_i * f_MBF_KK_i * f_MBF_Bs_i * X_JR_i * X_JB_i *
                         fcn_i;
        EvtComplex a_st_i = conj( a_i );
        EvtComplex a_f = f_PHSP_f * f_MBF_KK_f * f_MBF_Bs_f * X_JR_f * X_JB_f *
                         fcn_f;
        EvtComplex a_st_f = conj( a_f );

        integral += 0.5 * bin_width * ( a_i * a_st_i + a_f * a_st_f );
    }

    return sqrt( abs2( integral ) );
}

// Blatt-Weisskopf barrier factors
double EvtBsMuMuKK::X_J( const int J, const double q, const int isB ) const
{
    double r_BW = 1.0;

    if ( isB == 0 ) {
        r_BW = 1.5;
    } else if ( isB == 1 ) {
        r_BW = 5.0;
    }

    double zsq = pow( r_BW * q, 2 );

    double X_J( 1.0 );

    if ( J == 1 ) {
        X_J = sqrt( 1.0 / ( 1.0 + zsq ) );
    } else if ( J == 2 ) {
        X_J = sqrt( 1.0 / ( zsq * zsq + 3.0 * zsq + 9.0 ) );
    }

    return X_J;
}

// EvtGen d matrix: Input is 2J instead of J etc
double EvtBsMuMuKK::Wignerd( const int J, const int l, const int alpha,
                             const double theta ) const
{
    return EvtdFunction::d( 2 * J, 2 * l, 2 * alpha, theta );
}

// J spin of KK, l spin proj of J/psi, alpha dimuon spin
EvtComplex EvtBsMuMuKK::AngularDist( const int J, const int l, const int alpha,
                                     const double cK, const double cL,
                                     const double chi ) const
{
    double thetaL = acos( cL );
    double thetaK = acos( cK );

    EvtComplex out = 0.5 * sqrt( ( 2 * J + 1 ) / pi ) *
                     exp( EvtComplex( 0, -l * chi ) );

    out *= Wignerd( 1, l, alpha, thetaL ) * Wignerd( J, -l, 0, thetaK );

    return out;
}

// Time-dependent amplitude calculation
EvtComplex EvtBsMuMuKK::AmpTime( const int q, const EvtComplex& gplus,
                                 const EvtComplex& gminus, const double delta,
                                 const double lambda_abs, const double Amp,
                                 const double phis, const int eta ) const
{
    EvtComplex amp_time = Amp * EvtComplex( cos( -delta ), sin( -delta ) );
    double qphis = q * phis;
    amp_time *= ( gplus + eta * pow( lambda_abs, -1.0 * q ) *
                              EvtComplex( cos( qphis ), sin( qphis ) ) * gminus );

    if ( q == 1 ) {
        amp_time *= eta;
    }

    return amp_time;
}
