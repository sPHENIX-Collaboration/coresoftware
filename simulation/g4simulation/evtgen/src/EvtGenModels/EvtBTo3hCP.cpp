
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

#include "EvtGenModels/EvtBTo3hCP.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <cmath>
#include <stdlib.h>
#include <string>

#define square( x ) ( ( x ) * ( x ) )

/*
 * MK - 25/Aug/2016
 * The code bellow is not necessarilly well writen as it is 1-to-1 rewrite
 * from FORTRAN code (which is spread over 4 different source codes in not
 * fully obvious way). Once it is rewriten and giving correct results, I
 * will think how to give it more proper structure (mainly by checking what
 * is duplicated and how to simplify it).
 */

void EvtBTo3hCP::setConstants( double balpha, double bbeta )
{
    alphaCP = balpha;
    double calpha = cos( alphaCP );
    double salpha = sin( alphaCP );
    betaCP = bbeta;
    double cbeta = cos( betaCP );
    double sbeta = sin( betaCP );

    MA2 = square( M_B ) + square( M_pip ) + square( M_pi0 ) + square( M_pi0 );
    MB2 = square( M_B ) + square( M_pip ) + square( M_pim ) + square( M_pi0 );
    MC2 = square( M_B ) + square( M_Kp ) + square( M_pim ) + square( M_pi0 );

    double StrongPhase = 0;
    EvtComplex StrongExp( cos( StrongPhase ), sin( StrongPhase ) );

    EvtComplex Mat_Tp0( calpha, -salpha );
    Mat_Tp0 *= 1.09;
    EvtComplex Mat_Tm0( calpha, -salpha );
    Mat_Tm0 *= 1.09;
    EvtComplex Mat_T0p( calpha, -salpha );
    Mat_T0p *= 0.66;
    EvtComplex Mat_T0m( calpha, -salpha );
    Mat_T0m *= 0.66;
    EvtComplex Mat_Tpm( calpha, -salpha );
    Mat_Tpm *= 1.00;
    EvtComplex Mat_Tmp( calpha, -salpha );
    Mat_Tmp *= 0.47;

    EvtComplex Mat_Ppm( cos( -0.5 ), sin( -0.5 ) );
    Mat_Ppm *= -0.2;
    EvtComplex Mat_Pmp( cos( 2. ), sin( 2. ) );
    Mat_Pmp *= 0.15;

    EvtComplex Mat_P1 = 0.5 * ( Mat_Ppm - Mat_Pmp );
    EvtComplex Mat_P0 = 0.5 * ( Mat_Ppm + Mat_Pmp );
    EvtComplex Nat_Tp0( calpha, salpha );
    Nat_Tp0 *= 1.09;
    EvtComplex Nat_Tm0( calpha, salpha );
    Nat_Tm0 *= 1.09;
    EvtComplex Nat_T0p( calpha, salpha );
    Nat_T0p *= 0.66;
    EvtComplex Nat_T0m( calpha, salpha );
    Nat_T0m *= 0.66;
    EvtComplex Nat_Tpm( calpha, salpha );
    Nat_Tpm *= 1.00;
    EvtComplex Nat_Tmp( calpha, salpha );
    Nat_Tmp *= 0.47;
    EvtComplex Nat_P1 = Mat_P1;
    EvtComplex Nat_P0 = Mat_P0;

    Mat_Tpm = StrongExp * Mat_Tpm;
    Nat_Tpm = StrongExp * Nat_Tpm;

    Mat_S1 = Mat_Tp0 + 2. * Mat_P1;
    Mat_S2 = Mat_T0p - 2. * Mat_P1;
    Mat_S3 = Mat_Tpm + Mat_P1 + Mat_P0;
    Mat_S4 = Mat_Tmp - Mat_P1 + Mat_P0;
    Mat_S5 = -Mat_Tpm - Mat_Tmp + Mat_Tp0 + Mat_T0p - 2. * Mat_P0;

    Nat_S1 = Nat_Tp0 + 2. * Nat_P1;
    Nat_S2 = Nat_T0p - 2. * Nat_P1;
    Nat_S3 = Nat_Tpm + Nat_P1 + Nat_P0;
    Nat_S4 = Nat_Tmp - Nat_P1 + Nat_P0;
    Nat_S5 = -Nat_Tpm - Nat_Tmp + Nat_Tp0 + Nat_T0p - 2. * Nat_P0;

    // B0    -->-- K*+ pi- Amplitudes (Trees + Penguins)
    MatKstarp = EvtComplex( calpha, -salpha ) * EvtComplex( 0.220, 0. ) +
                EvtComplex( cbeta, sbeta ) * EvtComplex( -1.200, 0. );
    // B0    -->-- K*0 pi0 Amplitudes (Trees + Penguins)
    MatKstar0 = EvtComplex( calpha, -salpha ) * EvtComplex( 0.015, 0. ) +
                EvtComplex( cbeta, sbeta ) * EvtComplex( 0.850, 0. );
    // B0    -->-- K+ rho- Amplitudes (Trees + Penguins)
    MatKrho = EvtComplex( calpha, -salpha ) * EvtComplex( 0.130, 0. ) +
              EvtComplex( cbeta, sbeta ) * EvtComplex( 0.160, 0. );
    // B0bar -->-- K*+ pi- Amplitudes (Trees + Penguins)
    NatKstarp = EvtComplex( 0., 0. );
    // B0bar -->-- K*0 pi0 Amplitudes (Trees + Penguins)
    NatKstar0 = EvtComplex( 0., 0. );
    // B0bar -->-- K+ rho- Amplitudes (Trees + Penguins)
    NatKrho = EvtComplex( 0., 0. );

}

void EvtBTo3hCP::Evt3pi( double alpha, int iset, EvtVector4R& p_pi_plus,
                         EvtVector4R& p_pi_minus, EvtVector4R& p_gamma_1,
                         EvtVector4R& p_gamma_2, double& Real_B0,
                         double& Imag_B0, double& Real_B0bar, double& Imag_B0bar )
{
    EvtVector4R p_p2;
    double AB0, AB0bar, Ainter, R1, R2;
    double factor;
    int ierr = 0;

    // Ghm : beta is not needed for this generation - put a default value
    setConstants( alpha, 0.362 );

    if ( iset == 0 ) {
        p_pi_plus.set( M_pip, 0, 0, 0 );
        p_p2.set( M_pi0, 0, 0, 0 );
        p_pi_minus.set( M_pim, 0, 0, 0 );

        do {
            firstStep( p_pi_plus, p_p2, p_pi_minus, 1 );
            ierr = compute3pi( p_pi_plus, p_p2, p_pi_minus, Real_B0, Imag_B0,
                               Real_B0bar, Imag_B0bar, iset );
        } while ( ierr != 0 );
    } else if ( iset < 0 ) {
        p_p2 = p_gamma_1 + p_gamma_2;
        ierr = compute3pi( p_pi_plus, p_p2, p_pi_minus, Real_B0, Imag_B0,
                           Real_B0bar, Imag_B0bar, iset );
        if ( ierr != 0 ) {
            std::cout << "Provided kinematics is not physical\n";
            std::cout << "Program will stop\n";
            exit( 1 );
        }
    } else    // iset > 0
    {
        factor_max = 0;

        int endLoop = iset;
        for ( int i = 0; i < endLoop; ++i ) {
            p_pi_plus.set( M_pip, 0, 0, 0 );
            p_p2.set( M_pi0, 0, 0, 0 );
            p_pi_minus.set( M_pim, 0, 0, 0 );

            firstStep( p_pi_plus, p_p2, p_pi_minus, 1 );
            ierr = compute3pi( p_pi_plus, p_p2, p_pi_minus, Real_B0, Imag_B0,
                               Real_B0bar, Imag_B0bar, iset );
            if ( ierr != 0 ) {
                continue;
            }
            AB0 = square( Real_B0 ) + square( Imag_B0 );
            AB0bar = square( Real_B0bar ) + square( Imag_B0bar );
            Ainter = Real_B0 * Imag_B0bar - Imag_B0 * Real_B0bar;
            R1 = ( AB0 - AB0bar ) / ( AB0 + AB0bar );
            R2 = ( 2.0 * Ainter ) / ( AB0 + AB0bar );
            factor = ( 1.0 + sqrt( square( R1 ) + square( R2 ) ) ) *
                     ( AB0 + AB0bar ) / 2.0;
            if ( factor > factor_max )
                factor_max = factor;
        }
        factor_max = 1.0 / std::sqrt( factor_max );
    }

    Real_B0 *= factor_max;
    Imag_B0 *= factor_max;
    Real_B0bar *= factor_max;
    Imag_B0bar *= factor_max;

    if ( iset < 0 ) {
        return;
    }

    rotation( p_pi_plus, 1 );
    rotation( p_p2, 0 );
    rotation( p_pi_minus, 0 );

    gammaGamma( p_p2, p_gamma_1, p_gamma_2 );
}

void EvtBTo3hCP::Evt3piMPP( double alpha, int iset, EvtVector4R& p_p1,
                            EvtVector4R& p_p2, EvtVector4R& p_p3,
                            double& Real_B0, double& Imag_B0,
                            double& Real_B0bar, double& Imag_B0bar )
{
    double ABp, ABm;
    int ierr = 0;

    // Ghm : beta is not needed for this generation - put a default value
    setConstants( alpha, 0.362 );

    if ( iset == 0 ) {
        p_p1.set( M_pim, 0, 0, 0 );
        p_p2.set( M_pip, 0, 0, 0 );
        p_p3.set( M_pip, 0, 0, 0 );

        do {
            firstStep( p_p1, p_p2, p_p3, 2 );
            ierr = compute3piMPP( p_p1, p_p2, p_p3, Real_B0, Imag_B0,
                                  Real_B0bar, Imag_B0bar, iset );
        } while ( ierr != 0 );
    } else if ( iset < 0 ) {
        ierr = compute3piMPP( p_p1, p_p2, p_p3, Real_B0, Imag_B0, Real_B0bar,
                              Imag_B0bar, iset );
        if ( ierr != 0 ) {
            std::cout << "Provided kinematics is not physical\n";
            std::cout << "Program will stop\n";
            exit( 1 );
        }
    } else    // iset > 0
    {
        factor_max = 0;

        int endLoop = iset;
        for ( int i = 0; i < endLoop; ++i ) {
            p_p1.set( M_pim, 0, 0, 0 );
            p_p2.set( M_pip, 0, 0, 0 );
            p_p3.set( M_pip, 0, 0, 0 );

            firstStep( p_p1, p_p2, p_p3, 2 );
            ierr = compute3piMPP( p_p1, p_p2, p_p3, Real_B0, Imag_B0,
                                  Real_B0bar, Imag_B0bar, iset );
            if ( ierr != 0 ) {
                continue;
            }
            ABp = square( Real_B0 ) + square( Imag_B0 );
            ABm = square( Real_B0bar ) + square( Imag_B0bar );
            if ( ABp > factor_max )
                factor_max = ABp;
            if ( ABm > factor_max )
                factor_max = ABm;
        }
        factor_max = 1.0 / std::sqrt( factor_max );
    }

    Real_B0 *= factor_max;
    Imag_B0 *= factor_max;
    Real_B0bar *= factor_max;
    Imag_B0bar *= factor_max;

    if ( iset < 0 ) {
        return;
    }

    rotation( p_p1, 1 );
    rotation( p_p2, 0 );
    rotation( p_p3, 0 );
}

void EvtBTo3hCP::Evt3piP00( double alpha, int iset, EvtVector4R& p_p1,
                            EvtVector4R& p_p1_gamma1, EvtVector4R& p_p1_gamma2,
                            EvtVector4R& p_p2_gamma1, EvtVector4R& p_p2_gamma2,
                            double& Real_B0, double& Imag_B0,
                            double& Real_B0bar, double& Imag_B0bar )
{
    double ABp, ABm;
    EvtVector4R p_p2, p_p3;
    int ierr = 0;

    // Ghm : beta is not needed for this generation - put a default value
    setConstants( alpha, 0.362 );

    if ( iset == 0 ) {
        p_p1.set( M_pip, 0, 0, 0 );
        p_p2.set( M_pi0, 0, 0, 0 );
        p_p3.set( M_pi0, 0, 0, 0 );

        do {
            firstStep( p_p1, p_p2, p_p3, 3 );
            ierr = compute3piP00( p_p1, p_p2, p_p3, Real_B0, Imag_B0,
                                  Real_B0bar, Imag_B0bar, iset );
        } while ( ierr != 0 );
    } else if ( iset < 0 ) {
        p_p2 = p_p1_gamma1 + p_p1_gamma2;
        p_p3 = p_p2_gamma1 + p_p2_gamma2;
        ierr = compute3piP00( p_p1, p_p2, p_p3, Real_B0, Imag_B0, Real_B0bar,
                              Imag_B0bar, iset );
        if ( ierr != 0 ) {
            std::cout << "Provided kinematics is not physical\n";
            std::cout << "Program will stop\n";
            exit( 1 );
        }
    } else    // iset > 0
    {
        factor_max = 0;

        int endLoop = iset;
        for ( int i = 0; i < endLoop; ++i ) {
            p_p1.set( M_pip, 0, 0, 0 );
            p_p2.set( M_pi0, 0, 0, 0 );
            p_p3.set( M_pi0, 0, 0, 0 );

            firstStep( p_p1, p_p2, p_p3, 3 );
            ierr = compute3piP00( p_p1, p_p2, p_p3, Real_B0, Imag_B0,
                                  Real_B0bar, Imag_B0bar, iset );
            if ( ierr != 0 ) {
                continue;
            }
            ABp = square( Real_B0 ) + square( Imag_B0 );
            ABm = square( Real_B0bar ) + square( Imag_B0bar );
            if ( ABp > factor_max )
                factor_max = ABp;
            if ( ABm > factor_max )
                factor_max = ABm;
        }
        factor_max = 1.0 / std::sqrt( factor_max );
    }

    Real_B0 *= factor_max;
    Imag_B0 *= factor_max;
    Real_B0bar *= factor_max;
    Imag_B0bar *= factor_max;

    if ( iset < 0 ) {
        return;
    }

    rotation( p_p1, 1 );
    rotation( p_p2, 0 );
    rotation( p_p3, 0 );

    gammaGamma( p_p2, p_p1_gamma1, p_p1_gamma2 );
    gammaGamma( p_p3, p_p2_gamma1, p_p2_gamma2 );
}

void EvtBTo3hCP::EvtKpipi( double alpha, double beta, int iset,
                           EvtVector4R& p_K_plus, EvtVector4R& p_pi_minus,
                           EvtVector4R& p_gamma_1, EvtVector4R& p_gamma_2,
                           double& Real_B0, double& Imag_B0, double& Real_B0bar,
                           double& Imag_B0bar )
{
    EvtVector4R p_p3;
    double ABp, ABm;
    int ierr = 0;

    setConstants( alpha, beta );

    if ( iset == 0 ) {
        p_K_plus.set( M_Kp, 0, 0, 0 );
        p_pi_minus.set( M_pim, 0, 0, 0 );
        p_p3.set( M_pi0, 0, 0, 0 );

        do {
            firstStep( p_K_plus, p_pi_minus, p_p3, 0 );
            ierr = computeKpipi( p_K_plus, p_pi_minus, p_p3, Real_B0, Imag_B0,
                                 Real_B0bar, Imag_B0bar, iset );
        } while ( ierr != 0 );
    } else if ( iset < 0 ) {
        p_p3 = p_gamma_1 + p_gamma_2;
        ierr = computeKpipi( p_K_plus, p_pi_minus, p_p3, Real_B0, Imag_B0,
                             Real_B0bar, Imag_B0bar, iset );
        if ( ierr != 0 ) {
            std::cout << "Provided kinematics is not physical\n";
            std::cout << "Program will stop\n";
            exit( 1 );
        }
    } else    // iset > 0
    {
        factor_max = 0;

        int endLoop = iset;
        for ( int i = 0; i < endLoop; ++i ) {
            p_K_plus.set( M_Kp, 0, 0, 0 );
            p_pi_minus.set( M_pim, 0, 0, 0 );
            p_p3.set( M_pi0, 0, 0, 0 );
            firstStep( p_K_plus, p_pi_minus, p_p3, 0 );
            ierr = computeKpipi( p_K_plus, p_pi_minus, p_p3, Real_B0, Imag_B0,
                                 Real_B0bar, Imag_B0bar, iset );
            if ( ierr != 0 ) {
                continue;
            }
            ABp = square( Real_B0 ) + square( Imag_B0 );
            ABm = square( Real_B0bar ) + square( Imag_B0bar );

            if ( ABp > factor_max ) {
                factor_max = ABp;
            }
            if ( ABm > factor_max ) {
                factor_max = ABm;
            }
        }
        factor_max = 1.0 / std::sqrt( factor_max );
    }

    Real_B0 *= factor_max;
    Imag_B0 *= factor_max;
    Real_B0bar *= factor_max;
    Imag_B0bar *= factor_max;

    if ( iset < 0 ) {
        return;
    }

    rotation( p_K_plus, 1 );
    rotation( p_pi_minus, 0 );
    rotation( p_p3, 0 );

    gammaGamma( p_p3, p_gamma_1, p_gamma_2 );
}

void EvtBTo3hCP::firstStep( EvtVector4R& p1, EvtVector4R& p2, EvtVector4R& p3,
                            int mode )
{
    const double m1sq = p1.mass2();
    const double m2sq = p2.mass2();
    const double m3sq = p3.mass2();
    double min_m12, min_m13, min_m23;

    double max_m12 = square( M_B );
    double max_m13 = square( M_B );
    double max_m23 = square( M_B );

    if ( mode == 0 ) {
        min_m12 = m1sq + m2sq + 2 * sqrt( m1sq * m2sq );
        min_m13 = m1sq + m3sq + 2 * sqrt( m1sq * m3sq );
        min_m23 = m2sq + m3sq + 2 * sqrt( m2sq * m3sq );
    } else {
        min_m12 = m1sq + m2sq;
        min_m13 = m1sq + m3sq;
        min_m23 = m2sq + m3sq;
    }

    bool eventOK;
    double m13, m12, m23;
    double E1;
    double E2;
    double E3;
    double p1mom;
    double p2mom;
    double p3mom;
    double cost13;
    double cost12;
    double cost23;
    eventOK = false;

    do {
        switch ( mode ) {
            case 0:
                generateSqMasses_Kpipi( m12, m13, m23, MC2, m1sq, m2sq, m3sq );
                break;
            case 1:
                generateSqMasses_3pi( m12, m13, m23, MB2, m1sq, m2sq, m3sq );
                break;
            case 2:
                generateSqMasses_3piMPP( m12, m13, m23, MB2, m1sq, m2sq, m3sq );
                break;
            case 3:
                generateSqMasses_3piP00( m12, m13, m23, MA2, m1sq, m2sq, m3sq );
                break;
            default:
                break;
        }
        // Check whether event is physical
        if ( ( m23 < min_m23 ) || ( m23 > max_m23 ) )
            continue;
        if ( ( m13 < min_m13 ) || ( m13 > max_m13 ) )
            continue;
        if ( ( m12 < min_m12 ) || ( m12 > max_m12 ) )
            continue;

        // Now check the cosines of the angles
        E1 = ( square( M_B ) + m1sq - m23 ) / ( 2. * M_B );
        E2 = ( square( M_B ) + m2sq - m13 ) / ( 2. * M_B );
        E3 = ( square( M_B ) + m3sq - m12 ) / ( 2. * M_B );
        p1mom = square( E1 ) - m1sq;
        p2mom = square( E2 ) - m2sq;
        p3mom = square( E3 ) - m3sq;
        if ( p1mom < 0 || p2mom < 0 || p3mom < 0 ) {
            //      std::cout<<"Momenta magnitude negative\n";
            continue;
        }
        p1mom = sqrt( p1mom );
        p2mom = sqrt( p2mom );
        p3mom = sqrt( p3mom );

        cost13 = ( 2. * E1 * E3 + m1sq + m3sq - m13 ) / ( 2. * p1mom * p3mom );
        cost12 = ( 2. * E1 * E2 + m1sq + m2sq - m12 ) / ( 2. * p1mom * p2mom );
        cost23 = ( 2. * E2 * E3 + m2sq + m3sq - m23 ) / ( 2. * p2mom * p3mom );
        if ( cost13 < -1. || cost13 > 1. || cost12 < -1. || cost12 > 1. ||
             cost23 < -1. || cost23 > 1. ) {
            continue;
        }
        eventOK = true;
    } while ( eventOK == false );

    // Now is time to fill 4-vectors
    p3.set( E3, 0, 0, p3mom );
    p1.set( E1, p1mom * sqrt( 1 - square( cost13 ) ), 0, p1mom * cost13 );
    p2.set( 0, E2 );
    for ( int i = 1; i < 4; ++i ) {
        p2.set( i, -p1.get( i ) - p3.get( i ) );
    }
    if ( p1.get( 0 ) < p1.d3mag() ) {
        std::cout << "Unphysical p1 generated: " << p1 << std::endl;
    }
    if ( p2.get( 0 ) < p2.d3mag() ) {
        std::cout << "Unphysical p2 generated: " << p2 << std::endl;
    }
    if ( p3.get( 0 ) < p3.d3mag() ) {
        std::cout << "Unphysical p3 generated: " << p3 << std::endl;
    }
    double testMB2 = MB2;
    switch ( mode ) {
        case 0:
            testMB2 = MC2;
            break;
        case 1:
        case 2:
            testMB2 = MB2;
            break;
        case 3:
            testMB2 = MA2;
            break;
    }

    if ( fabs( m12 + m13 + m23 - testMB2 ) > 1e-4 ) {
        std::cout << "Unphysical event generated: " << m12 << " " << m13 << " "
                  << m23 << std::endl;
    }
}

void EvtBTo3hCP::generateSqMasses_Kpipi( double& m12, double& m13, double& m23,
                                         double MB2, double m1sq, double m2sq,
                                         double m3sq )
{
    /*
  C There is two ways of generating the events:
  C The first one used a pole-compensation method to generate the
  C events efficiently taking into account the poles due to the
  C Breit-Wigners of the rho s. It is activated by setting
  C Phase_Space to .false.
  C The second one generates events according to phase space. It is
  C inneficient but allows the exploration of the full Dalitz plot
  C in an uniform way. It was found to be usefull fopr some peculiar
  C applications. It is activated by setting
  C Phase_space to .true.
  C Note that in that case, the generation is no longer correct.
  */
    static bool phaseSpace = false;

    double max_m12 = square( M_B );
    double min_m12 = m1sq + m2sq + 2 * sqrt( m1sq * m2sq );

    double max_m13 = square( M_B );
    double min_m13 = m1sq + m3sq + 2 * sqrt( m1sq * m3sq );

    double max_m23 = square( M_B );
    double min_m23 = m2sq + m3sq + 2 * sqrt( m2sq * m3sq );

    double z = 3. * EvtRandom::Flat();
    if ( z < 1. )    // K*+
    {
        if ( phaseSpace ) {
            m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        } else {
            double y = EvtRandom::Flat() * pi - pi / 2;
            double x = std::tan( y );
            double mass = x * Gam_Kstarp / 2. + Mass_Kstarp;
            m13 = square( mass );
        }
        m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        m23 = MB2 - m12 - m13;
    } else if ( z < 2. )    // K*0
    {
        if ( phaseSpace ) {
            m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        } else {
            double y = EvtRandom::Flat() * pi - pi / 2;
            double x = std::tan( y );
            double mass = x * Gam_Kstar0 / 2. + Mass_Kstar0;
            m12 = square( mass );
        }
        m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        m23 = MB2 - m12 - m13;
    } else    // rho-
    {
        if ( phaseSpace ) {
            m23 = EvtRandom::Flat() * ( max_m23 - min_m23 ) + min_m23;
        } else {
            double y = EvtRandom::Flat() * pi - pi / 2;
            double x = std::tan( y );
            double mass = x * Gam_rho / 2. + Mass_rho;
            m23 = square( mass );
        }
        m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        m12 = MB2 - m23 - m13;
    }
}

void EvtBTo3hCP::generateSqMasses_3pi( double& m12, double& m13, double& m23,
                                       double MB2, double m1sq, double m2sq,
                                       double m3sq )
{
    /*
  C There is two ways of generating the events:
  C The first one used a pole-compensation method to generate the
  C events efficiently taking into account the poles due to the
  C Breit-Wigners of the rho s. It is activated by setting
  C Phase_Space to .false.
  C The second one generates events according to phase space. It is
  C inneficient but allows the exploration of the full Dalitz plot
  C in an uniform way. It was found to be usefull fopr some peculiar
  C applications. It is activated by setting
  C Phase_space to .true.
  C Note that in that case, the generation is no longer correct.
  */
    static bool phaseSpace = false;

    double max_m12 = square( M_B );
    double min_m12 = m1sq + m2sq;

    double max_m13 = square( M_B );
    double min_m13 = m1sq + m3sq;

    double max_m23 = square( M_B );
    double min_m23 = m2sq + m3sq;
    double mass = 0;

    if ( !phaseSpace ) {
        double y = EvtRandom::Flat() * pi - pi / 2;
        double x = std::tan( y );
        mass = x * Gam_rho / 2. + Mass_rho;
    }

    double z = 3. * EvtRandom::Flat();
    if ( z < 1. ) {
        if ( phaseSpace ) {
            m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        } else {
            m12 = square( mass );
        }
        m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        m23 = MB2 - m12 - m13;
    } else if ( z < 2. ) {
        if ( phaseSpace ) {
            m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        } else {
            m13 = square( mass );
        }
        m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        m23 = MB2 - m12 - m13;
    } else {
        if ( phaseSpace ) {
            m23 = EvtRandom::Flat() * ( max_m23 - min_m23 ) + min_m23;
        } else {
            m23 = square( mass );
        }
        m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        m13 = MB2 - m12 - m23;
    }
}

void EvtBTo3hCP::generateSqMasses_3piMPP( double& m12, double& m13, double& m23,
                                          double MB2, double m1sq, double m2sq,
                                          double m3sq )
{
    /*
  C There is two ways of generating the events:
  C The first one used a pole-compensation method to generate the
  C events efficiently taking into account the poles due to the
  C Breit-Wigners of the rho s. It is activated by setting
  C Phase_Space to .false.
  C The second one generates events according to phase space. It is
  C inneficient but allows the exploration of the full Dalitz plot
  C in an uniform way. It was found to be usefull fopr some peculiar
  C applications. It is activated by setting
  C Phase_space to .true.
  C Note that in that case, the generation is no longer correct.
  */
    static bool phaseSpace = false;

    double max_m12 = square( M_B );
    double min_m12 = m1sq + m2sq;

    double max_m13 = square( M_B );
    double min_m13 = m1sq + m3sq;

    double mass = 0;

    if ( !phaseSpace ) {
        double y = EvtRandom::Flat() * pi - pi / 2;
        double x = std::tan( y );
        mass = x * Gam_rho / 2. + Mass_rho;
    }

    double z = EvtRandom::Flat();
    if ( z < 0.5 ) {
        if ( phaseSpace ) {
            m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        } else {
            m12 = square( mass );
        }
        m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        m23 = MB2 - m12 - m13;
    } else {
        if ( phaseSpace ) {
            m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        } else {
            m13 = square( mass );
        }
        m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        m23 = MB2 - m12 - m13;
    }
}

void EvtBTo3hCP::generateSqMasses_3piP00( double& m12, double& m13, double& m23,
                                          double MB2, double m1sq, double m2sq,
                                          double m3sq )
{
    /*
  C There is two ways of generating the events:
  C The first one used a pole-compensation method to generate the
  C events efficiently taking into account the poles due to the
  C Breit-Wigners of the rho s. It is activated by setting
  C Phase_Space to .false.
  C The second one generates events according to phase space. It is
  C inneficient but allows the exploration of the full Dalitz plot
  C in an uniform way. It was found to be usefull fopr some peculiar
  C applications. It is activated by setting
  C Phase_space to .true.
  C Note that in that case, the generation is no longer correct.
  */
    static bool phaseSpace = false;

    double max_m12 = square( M_B );
    double min_m12 = m1sq + m2sq;

    double max_m13 = square( M_B );
    double min_m13 = m1sq + m3sq;

    double mass = 0;

    if ( !phaseSpace ) {
        double y = EvtRandom::Flat() * pi - pi / 2;
        double x = std::tan( y );
        mass = x * Gam_rho / 2. + Mass_rho;
    }

    double z = EvtRandom::Flat();
    if ( z < 0.5 ) {
        if ( phaseSpace ) {
            m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        } else {
            m12 = square( mass );
        }
        m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        m23 = MB2 - m12 - m13;
    } else {
        if ( phaseSpace ) {
            m13 = EvtRandom::Flat() * ( max_m13 - min_m13 ) + min_m13;
        } else {
            m13 = square( mass );
        }
        m12 = EvtRandom::Flat() * ( max_m12 - min_m12 ) + min_m12;
        m23 = MB2 - m12 - m13;
    }
}
int EvtBTo3hCP::compute3pi( EvtVector4R& p1, EvtVector4R& p2, EvtVector4R& p3,
                            double& real_B0, double& imag_B0,
                            double& real_B0bar, double& imag_B0bar, int iset )
{
    int ierr = 0;

    double m12 = ( p1 + p2 ).mass();
    double m13 = ( p1 + p3 ).mass();
    double m23 = ( p2 + p3 ).mass();

    double W12 = 1. /
                 ( ( square( Mass_rho - m12 ) + square( Gam_rho / 2. ) ) * m12 );
    double W13 = 1. /
                 ( ( square( Mass_rho - m13 ) + square( Gam_rho / 2. ) ) * m13 );
    double W23 = 1. /
                 ( ( square( Mass_rho - m23 ) + square( Gam_rho / 2. ) ) * m23 );

    double Wtot = 1.;
    if ( iset >= 0 ) {
        Wtot = 1. / sqrt( W12 + W13 + W23 );
    }

    EvtComplex Mat_rhop = BreitWigner( p1, p2, p3, ierr );
    EvtComplex Mat_rhom = BreitWigner( p2, p3, p1, ierr );
    EvtComplex Mat_rho0 = BreitWigner( p1, p3, p2, ierr );

    EvtComplex Mat_1 = Mat_S3 * Mat_rhop;
    EvtComplex Mat_2 = Mat_S4 * Mat_rhom;
    EvtComplex Mat_3 = Mat_S5 * Mat_rho0 * 0.5;

    EvtComplex MatBp = ( Mat_1 + Mat_2 + Mat_3 ) * Wtot;

    Mat_1 = Nat_S3 * Mat_rhom;
    Mat_2 = Nat_S4 * Mat_rhop;
    Mat_3 = Nat_S5 * Mat_rho0 * 0.5;

    EvtComplex MatBm = ( Mat_1 + Mat_2 + Mat_3 ) * Wtot;

    real_B0 = real( MatBp );
    imag_B0 = imag( MatBp );

    real_B0bar = real( MatBm );
    imag_B0bar = imag( MatBm );

    return ierr;
}

int EvtBTo3hCP::compute3piMPP( EvtVector4R& p1, EvtVector4R& p2,
                               EvtVector4R& p3, double& real_B0, double& imag_B0,
                               double& real_B0bar, double& imag_B0bar, int iset )
{
    int ierr = 0;
    const double ASHQ = sqrt( 2. );
    double m12 = ( p1 + p2 ).mass();
    double m13 = ( p1 + p3 ).mass();

    double W12 = 1. /
                 ( ( square( Mass_rho - m12 ) + square( Gam_rho / 2. ) ) * m12 );
    double W13 = 1. /
                 ( ( square( Mass_rho - m13 ) + square( Gam_rho / 2. ) ) * m13 );

    double Wtot = 1.;
    if ( iset >= 0 ) {
        Wtot = 1. / sqrt( W12 + W13 );
    }

    EvtComplex Mat_rhop = BreitWigner( p1, p2, p3, ierr ) +
                          BreitWigner( p1, p3, p2, ierr );

    EvtComplex MatBp = Mat_S2 * Mat_rhop * Wtot * ASHQ;
    EvtComplex MatBm = Nat_S2 * Mat_rhop * Wtot * ASHQ;

    real_B0 = real( MatBp );
    imag_B0 = imag( MatBp );

    real_B0bar = real( MatBm );
    imag_B0bar = imag( MatBm );

    return ierr;
}

int EvtBTo3hCP::compute3piP00( EvtVector4R& p1, EvtVector4R& p2,
                               EvtVector4R& p3, double& real_B0, double& imag_B0,
                               double& real_B0bar, double& imag_B0bar, int iset )
{
    int ierr = 0;
    const double ASHQ = sqrt( 2. );
    double m12 = ( p1 + p2 ).mass();
    double m13 = ( p1 + p3 ).mass();

    double W12 = 1. /
                 ( ( square( Mass_rho - m12 ) + square( Gam_rho / 2. ) ) * m12 );
    double W13 = 1. /
                 ( ( square( Mass_rho - m13 ) + square( Gam_rho / 2. ) ) * m13 );

    double Wtot = 1.;
    if ( iset >= 0 ) {
        Wtot = 1. / sqrt( W12 + W13 );
    }

    EvtComplex Mat_rhop = BreitWigner( p1, p2, p3, ierr ) +
                          BreitWigner( p1, p3, p2, ierr );

    EvtComplex MatBp = Mat_S1 * Mat_rhop * Wtot * ASHQ;
    EvtComplex MatBm = Nat_S1 * Mat_rhop * Wtot * ASHQ;

    real_B0 = real( MatBp );
    imag_B0 = imag( MatBp );

    real_B0bar = real( MatBm );
    imag_B0bar = imag( MatBm );

    return ierr;
}

int EvtBTo3hCP::computeKpipi( EvtVector4R& p1, EvtVector4R& p2, EvtVector4R& p3,
                              double& real_B0, double& imag_B0,
                              double& real_B0bar, double& imag_B0bar, int iset )
{
    int ierr = 0;

    double m12 = ( p1 + p2 ).mass();
    double m13 = ( p1 + p3 ).mass();
    double m23 = ( p2 + p3 ).mass();

    double W12 = 1. /
                 ( ( square( Mass_Kstar0 - m12 ) + square( Gam_Kstar0 / 2. ) ) *
                   m12 );
    double W13 = 1. /
                 ( ( square( Mass_Kstarp - m13 ) + square( Gam_Kstarp / 2. ) ) *
                   m13 );
    double W23 = 1. /
                 ( ( square( Mass_rho - m23 ) + square( Gam_rho / 2. ) ) * m23 );

    double Wtot = 1.;
    if ( iset >= 0 ) {
        Wtot = 1. / sqrt( W12 + W13 + W23 );
    }

    EvtComplex BW13 = BreitWigner( p1, p3, p2, ierr, Mass_Kstarp, Gam_Kstarp );
    if ( ierr != 0 )
        return ierr;
    EvtComplex BW12 = BreitWigner( p1, p2, p3, ierr, Mass_Kstar0, Gam_Kstar0 );
    if ( ierr != 0 )
        return ierr;
    /*
   If the rho is to be treated on the same footing as K* ==> use the line below
    EvtComplex BW23=BreitWigner(p2, p3, p1, ierr, Mass_Rho, Gam_Rho);
  */
    EvtComplex BW23 = BreitWigner( p2, p3, p1, ierr );
    if ( ierr != 0 )
        return ierr;

    // Build up amplitudes
    EvtComplex MatB0 = MatKstarp * BW13 + MatKstar0 * BW12 + MatKrho * BW23;
    EvtComplex MatB0bar = NatKstarp * BW13 + NatKstar0 * BW12 + NatKrho * BW23;

    real_B0 = real( MatB0 ) * Wtot;
    imag_B0 = imag( MatB0 ) * Wtot;
    real_B0bar = real( MatB0bar ) * Wtot;
    imag_B0bar = imag( MatB0bar ) * Wtot;

    return ierr;
}

void EvtBTo3hCP::rotation( EvtVector4R& p, int newRot )
{
    if ( newRot ) {
        double phi2 = EvtRandom::Flat() * 2. * pi;
        double phi3 = EvtRandom::Flat() * 2. * pi;

        double c1 = 2. * EvtRandom::Flat() - 1.;
        double c2 = cos( phi2 );
        double c3 = cos( phi3 );

        double s1 = sqrt( 1. - square( c1 ) );
        double s2 = sin( phi2 );
        double s3 = sin( phi3 );

        rotMatrix[0][0] = c1;
        rotMatrix[0][1] = s1 * c3;
        rotMatrix[0][2] = s1 * s3;
        rotMatrix[1][0] = -s1 * c2;
        rotMatrix[1][1] = c1 * c2 * c3 - s2 * s3;
        rotMatrix[1][2] = c1 * c2 * s3 + s2 * c3;
        rotMatrix[2][0] = s1 * s2;
        rotMatrix[2][1] = -c1 * s2 * c3 - c2 * s3;
        rotMatrix[2][2] = -c1 * s2 * s3 + c2 * c3;
    }

    double mom[3];
    for ( int i = 1; i < 4; ++i ) {
        mom[i - 1] = p.get( i );
        p.set( i, 0 );
    }
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
            p.set( i + 1, p.get( i + 1 ) + rotMatrix[i][j] * mom[j] );
        }
    }
}

void EvtBTo3hCP::gammaGamma( EvtVector4R& p, EvtVector4R& pgamma1,
                             EvtVector4R& pgamma2 )
{
    double EGammaCmsPi0 = sqrt( p.mass2() ) / 2.;

    double cosThetaRot = EvtRandom::Flat() * 2. - 1.;
    double sinThetaRot = sqrt( 1. - square( cosThetaRot ) );
    double PhiRot = EvtRandom::Flat() * 2. * pi;

    pgamma1.set( 1, EGammaCmsPi0 * sinThetaRot * cos( PhiRot ) );
    pgamma1.set( 2, EGammaCmsPi0 * sinThetaRot * sin( PhiRot ) );
    pgamma1.set( 3, EGammaCmsPi0 * cosThetaRot );
    pgamma1.set( 0, EGammaCmsPi0 );

    for ( int i = 1; i < 4; ++i ) {
        pgamma2.set( i, -pgamma1.get( i ) );
    }
    pgamma2.set( 0, pgamma1.get( 0 ) );

    pgamma1.applyBoostTo( p );
    pgamma2.applyBoostTo( p );
}

EvtComplex EvtBTo3hCP::BreitWigner( EvtVector4R& p1, EvtVector4R& p2,
                                    EvtVector4R& p3, int& ierr, double Mass,
                                    double Width )
{
    bool pipiMode = true;

    if ( Mass > 1e-5 ) {
        pipiMode = false;
    }

    bool relatBW = true;
    bool aleph = true;
    EvtComplex result( 0, 0 );
    ierr = 0;

    double m12 = ( p1 + p2 ).mass();
    double e12 = ( p1 + p2 ).get( 0 );
    double argu = 1. - square( m12 ) / square( e12 );
    double beta = 0;
    if ( argu > 0 ) {
        beta = sqrt( argu );
    } else {
        std::cout << "Abnormal beta ! Argu  = " << argu << std::endl;
        argu = 0;
    }
    double gamma = e12 / m12;

    double m13sq = ( p1 + p3 ).mass2();
    double costet = ( 2. * p1.get( 0 ) * p3.get( 0 ) - m13sq + p1.mass2() +
                      p3.mass2() ) /
                    ( 2. * p1.d3mag() * p3.d3mag() );

    double p1z = p1.d3mag() * costet;
    double p1zcms12 = gamma * ( p1z + beta * p1.get( 0 ) );
    double e1cms12 = gamma * ( p1.get( 0 ) + beta * p1z );
    double p1cms12 = sqrt( square( e1cms12 ) - p1.mass2() );
    double coscms = p1zcms12 / p1cms12;

    if ( pipiMode ) {
        if ( aleph ) {
            double m12_2 = square( m12 );
            result = coscms * EvtCRhoF_W( m12_2 );
        } else {
            double factor = 2 * ( square( Mass_rho - m12 ) +
                                  square( 0.5 * Gam_rho ) );
            factor = coscms * Gam_rho / factor;
            double numReal = ( Mass_rho - m12 ) * factor;
            double numImg = 0.5 * Gam_rho * factor;
            result = EvtComplex( numReal, numImg );
        }
    } else {
        if ( relatBW ) {
            double Am2Min = p1.mass2() + p2.mass2() + 2 * p1.mass() * p2.mass();
            result = coscms *
                     EvtRBW( square( m12 ), square( Mass ), Width, Am2Min );
        } else {
            double factor = 2 * ( square( Mass - m12 ) + square( 0.5 * Width ) );
            factor = coscms * Width / factor;
            double numReal = ( Mass - m12 ) * factor;
            double numImg = 0.5 * Width * factor;
            result = EvtComplex( numReal, numImg );
        }
    }

    return result;
}

EvtComplex EvtBTo3hCP::EvtCRhoF_W( double s )
{
    const bool kuhn_santa = true;    // type of Breit-Wigner formula
                                     //  double lambda = 1.0;

    double AmRho, GamRho, AmRhoP, GamRhoP, beta, AmRhoPP, GamRhoPP, gamma;
    if ( kuhn_santa ) {
        //...rho(770)
        AmRho = 0.7734;
        GamRho = 0.1477;
        //...rho(1450)
        AmRhoP = 1.465;
        GamRhoP = 0.696;
        beta = -0.229;
        //...rho(1700)
        AmRhoPP = 1.760;
        GamRhoPP = 0.215;
        gamma = 0.075;
    } else {
        //...rho(770)
        AmRho = 0.7757;
        GamRho = 0.1508;
        //...rho(1450)
        AmRhoP = 1.448;
        GamRhoP = 0.503;
        beta = -0.161;
        //...rho(1700)
        AmRhoPP = 1.757;
        GamRhoPP = 0.237;
        gamma = 0.076;
    }

    EvtComplex result( 0, 0 );

    if ( kuhn_santa ) {
        result = ( EvtcBW_KS( s, square( AmRho ), GamRho ) +    //!...BW-rho( 770)
                   EvtcBW_KS( s, square( AmRhoP ), GamRhoP ) *
                       ( beta ) +    //!...BW-rho(1450)
                   EvtcBW_KS( s, square( AmRhoPP ), GamRhoPP ) *
                       ( gamma ) ) /    //!...BW-rho(1700)
                 ( 1. + beta + gamma );
    } else {
        result = ( EvtcBW_GS( s, square( AmRho ), GamRho ) +
                   EvtcBW_GS( s, square( AmRhoP ), GamRhoP ) * ( beta ) +
                   EvtcBW_GS( s, square( AmRhoPP ), GamRhoPP ) * ( gamma ) ) /
                 ( 1. + beta + gamma );
    }
    return result;
}

EvtComplex EvtBTo3hCP::EvtRBW( double s, double Am2, double Gam, double Am2Min )
{
    EvtComplex result( 0, 0 );

    if ( s < Am2Min ) {
        return result;
    }

    double tmp = ( ( s - Am2Min ) / ( Am2 - Am2Min ) );
    double G = Gam * ( Am2 / s ) * sqrt( square( tmp ) * tmp );
    double D = square( Am2 - s ) + s * square( G );
    double X = Am2 * ( Am2 - s );
    double Y = Am2 * sqrt( s ) * G;

    result = EvtComplex( X / D, Y / D );
    return result;
}

EvtComplex EvtBTo3hCP::EvtcBW_KS( double s, double Am2, double Gam )
{
    EvtComplex result( 0, 0 );
    const double AmPi2 = square( 0.13956995 );
    return EvtRBW( s, Am2, Gam, 4. * AmPi2 );
}

EvtComplex EvtBTo3hCP::EvtcBW_GS( double s, double Am2, double Gam )
{
    EvtComplex result( 0, 0 );
    const double AmPi2 = square( 0.13956995 );

    if ( s < 4. * AmPi2 ) {
        return result;
    }

    double tmp = ( ( s - 4. * AmPi2 ) / ( Am2 - 4. * AmPi2 ) );

    double G = Gam * ( Am2 / s ) * sqrt( square( tmp ) * tmp );
    double z1 = Am2 - s + Evtfs( s, Am2, Gam );
    double z2 = sqrt( s ) * G;
    double z3 = Am2 + d( Am2 ) * Gam * sqrt( Am2 );

    double X = z3 * z1;
    double Y = z3 * z2;
    double N = square( z1 ) + square( z2 );

    result = EvtComplex( X / N, Y / N );

    return result;
}

double EvtBTo3hCP::d( double AmRho2 )
{
    const double lpi = 3.141593;
    const double AmPi = 0.13956995;
    const double AmPi2 = square( AmPi );
    double AmRho = sqrt( AmRho2 );
    double k_AmRho2 = k( AmRho2 );
    double result = 3. / lpi * AmPi2 / square( k_AmRho2 ) *
                        log( ( AmRho + 2. * k_AmRho2 ) / ( 2. * AmPi ) ) +
                    AmRho / ( 2. * pi * k_AmRho2 ) -
                    AmPi2 * AmRho / ( pi * ( square( k_AmRho2 ) * k_AmRho2 ) );
    return result;
}

double EvtBTo3hCP::k( double s )
{
    const double AmPi2 = square( 0.13956995 );
    return 0.5 * sqrt( s - 4. * AmPi2 );
}

double EvtBTo3hCP::Evtfs( double s, double AmRho2, double GamRho )
{
    double k_s = k( s );
    double k_Am2 = k( AmRho2 );

    return GamRho * AmRho2 / ( square( k_Am2 ) * k_Am2 ) *
           ( square( k_s ) * ( h( s ) - h( AmRho2 ) ) +
             ( AmRho2 - s ) * square( k_Am2 ) * dh_ds( AmRho2 ) );
}

double EvtBTo3hCP::h( double s )
{
    const double pi = 3.141593;
    const double AmPi = 0.13956995;
    double sqrts = sqrt( s );
    double k_s = k( s );
    return 2. / pi * ( k_s / sqrts ) *
           log( ( sqrts + 2. * k_s ) / ( 2. * AmPi ) );
}

double EvtBTo3hCP::dh_ds( double s )
{
    const double pi = 3.141593;
    return h( s ) * ( 1. / ( 8. * square( k( s ) ) ) - 1. / ( 2 * s ) ) +
           1. / ( 2. * pi * s );
}
