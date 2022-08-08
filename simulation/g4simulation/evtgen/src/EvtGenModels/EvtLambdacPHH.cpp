
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

#include "EvtGenModels/EvtLambdacPHH.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <algorithm>
#include <cmath>
#include <utility>

EvtLambdacPHH::EvtLambdacPHH() :
    _d1( 0 ),
    _d2( 1 ),
    _d3( 3 ),
    _Nplusplus( 0.46 ),
    _Nplusminus( 1.0 ),
    _Nminusplus( 0.18 ),
    _Nminusminus( 0.94 ),
    _phiNplusplus( 3.48 ),
    _phiNplusminus( 0.00 ),
    _phiNminusplus( 0.75 ),
    _phiNminusminus( 1.13 ),
    _E1( 0.52 ),
    _phiE1( -1.01 ),
    _E2( 0.20 ),
    _phiE2( 2.35 ),
    _E3( 0.21 ),
    _phiE3( 3.46 ),
    _E4( 0.16 ),
    _phiE4( 5.29 ),
    _F1( 0.17 ),
    _phiF1( 4.98 ),
    _F2( 0.38 ),
    _phiF2( 4.88 ),
    _H1( 0.18 ),
    _phiH1( 5.93 ),
    _H2( 0.20 ),
    _phiH2( -0.06 ),
    _NRNorm( 1.0 ),
    _KstarNorm( 1.0 ),
    _DeltaNorm( 1.0 ),
    _LambdaNorm( 1.0 ),
    _KstarM( 0.890 ),
    _KstarW( 0.0498 ),
    _KstarR( 3.40 ),
    _DeltaM( 1.232 ),
    _DeltaW( 0.1120 ),
    _DeltaR( 5.22 ),
    _LambdaM( 1.520 ),
    _LambdaW( 0.0156 ),
    _LambdaR( 6.29 ),
    _Lambda_cR( 5.07 ),
    _zprime(),
    _p4_Lambda_c(),
    _zpMag( 0.0 ),
    _p4_Lambdac_Mag( 0.0 )
{
    // Fermilab E791 values from MINUIT fit arXiv:hep-ex/9912003v1
}

std::string EvtLambdacPHH::getName()
{
    return "LAMBDAC_PHH";
}

EvtDecayBase* EvtLambdacPHH::clone()
{
    return new EvtLambdacPHH;
}

bool compareId( const std::pair<EvtId, int>& left,
                const std::pair<EvtId, int>& right )
{
    // Compare id numbers to achieve the ordering K-, pi+ and p
    bool result( false );

    int leftPDGid = EvtPDL::getStdHep( left.first );
    int rightPDGid = EvtPDL::getStdHep( right.first );

    if ( leftPDGid < rightPDGid ) {
        result = true;
    }

    return result;
}

void EvtLambdacPHH::init()
{
    static EvtId KM = EvtPDL::getId( "K-" );
    static EvtId PIP = EvtPDL::getId( "pi+" );
    static EvtId LAMBDAC = EvtPDL::getId( "Lambda_c+" );
    static EvtId LAMBDACB = EvtPDL::getId( "anti-Lambda_c-" );
    static EvtId PROTON = EvtPDL::getId( "p+" );

    // check that there are 0 or 1 arguments and 3 daughters
    checkNArg( 0, 1 );
    checkNDaug( 3 );

    EvtId parnum = getParentId();
    checkSpinParent( EvtSpinType::DIRAC );
    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );

    std::vector<std::pair<EvtId, int>> daughters;
    if ( parnum == LAMBDAC ) {
        for ( int i = 0; i < 3; ++i ) {
            daughters.push_back( std::make_pair( getDaug( i ), i ) );
        }
    } else {
        for ( int i = 0; i < 3; ++i ) {
            daughters.push_back(
                std::make_pair( EvtPDL::chargeConj( getDaug( i ) ), i ) );
        }
    }

    // Sort daughters, they will end up in the order KM, PIP and PROTON
    std::sort( daughters.begin(), daughters.end(), compareId );

    if ( parnum == LAMBDAC || parnum == LAMBDACB ) {
        if ( daughters[0].first == KM && daughters[1].first == PIP &&
             daughters[2].first == PROTON ) {
            _d1 = daughters[0].second;
            _d2 = daughters[1].second;
            _d3 = daughters[2].second;
        }
    }

    // Find resonance dynamics normalisations
    calcNormalisations();

    // Print out expected fit fractions
    getFitFractions();
}

void EvtLambdacPHH::calcNormalisations()
{
    // Generate events uniform in the Lambda_c Dalitz plot and find the
    // normalisation integrals of the Breit-Wigner lineshapes

    // Lambda_c -> K- pi+ p
    int nDaug( 3 );
    EvtVector4R p4Daug[3];

    double mDaug[3] = {EvtPDL::getMeanMass( EvtPDL::getId( "K-" ) ),
                       EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ),
                       EvtPDL::getMeanMass( EvtPDL::getId( "p+" ) )};

    double norm[3] = {0.0, 0.0, 0.0};

    // sample size
    int N( 100000 );
    for ( int i = 0; i < N; i++ ) {
        double mParent = EvtPDL::getMass( EvtPDL::getId( "Lambda_c+" ) );
        EvtVector4R p0( mParent, 0.0, 0.0, 0.0 );

        // Generate uniform 4 momenta
        EvtGenKine::PhaseSpace( nDaug, mDaug, p4Daug, mParent );

        EvtResonance2 LambdacpKpi1( p0, p4Daug[0], p4Daug[1], 1.0, 0.0, _KstarW,
                                    _KstarM, 1, true, _KstarR,
                                    _Lambda_cR );    // K*0 -> K- and pi+; L = 1
        EvtResonance2 LambdacpKpi2( p0, p4Daug[2], p4Daug[1], 1.0, 0.0, _DeltaW,
                                    _DeltaM, 1, true, _DeltaR,
                                    _Lambda_cR );    // Delta++ -> p and pi+; L = 1
        EvtResonance2 LambdacpKpi3(
            p0, p4Daug[2], p4Daug[0], 1.0, 0.0, _LambdaW, _LambdaM, 2, true,
            _LambdaR, _Lambda_cR );    // Lambda(1520) -> K- and p; L = 2

        // Sum amplitude magnitude squared
        norm[0] += abs2( LambdacpKpi1.resAmpl() );
        norm[1] += abs2( LambdacpKpi2.resAmpl() );
        norm[2] += abs2( LambdacpKpi3.resAmpl() );
    }

    // Set normalisation lineshape multiplication factors
    double N0( N * 1.0 );

    // Scale NR to get sensible relative fit fractions
    _NRNorm = 1.0 / 3.0;
    // Set this using a decay file parameter if required
    if ( getNArg() > 1 ) {
        _NRNorm = getArg( 1 );
    }

    if ( norm[0] > 0.0 ) {
        _KstarNorm = sqrt( N0 / norm[0] );
    }
    if ( norm[1] > 0.0 ) {
        _DeltaNorm = sqrt( N0 / norm[1] );
    }
    if ( norm[2] > 0.0 ) {
        _LambdaNorm = sqrt( N0 / norm[2] );
    }
}

void EvtLambdacPHH::getFitFractions()
{
    // Generate events uniform in the Lambda_c Dalitz plot and find the
    // fit fractions for each resonance

    // Lambda_c -> K- pi+ p
    int nDaug( 3 );
    EvtVector4R p4Daug[3];

    double mDaug[3] = {EvtPDL::getMeanMass( EvtPDL::getId( "K-" ) ),
                       EvtPDL::getMeanMass( EvtPDL::getId( "pi+" ) ),
                       EvtPDL::getMeanMass( EvtPDL::getId( "p+" ) )};

    double FitFracTop[4] = {0.0, 0.0, 0.0, 0.0};
    double FitFracDenom = 0.0;

    // sample size
    int N( 100000 );
    for ( int i = 0; i < N; i++ ) {
        double mParent = EvtPDL::getMass( EvtPDL::getId( "Lambda_c+" ) );
        EvtVector4R p0( mParent, 0.0, 0.0, 0.0 );

        // Generate uniform 4 momenta
        EvtGenKine::PhaseSpace( nDaug, mDaug, p4Daug, mParent );

        EvtResonance2 LambdacpKpi0( p0, p4Daug[0], p4Daug[1], 1.0, 0.0, 0.0, 0.0,
                                    0, true, 0.0, 0.0 );    // Non resonant (NR)
        EvtResonance2 LambdacpKpi1( p0, p4Daug[0], p4Daug[1], 1.0, 0.0, _KstarW,
                                    _KstarM, 1, true, _KstarR,
                                    _Lambda_cR );    // K*0 -> K- and pi+; L = 1
        EvtResonance2 LambdacpKpi2( p0, p4Daug[2], p4Daug[1], 1.0, 0.0, _DeltaW,
                                    _DeltaM, 1, true, _DeltaR,
                                    _Lambda_cR );    // Delta++ -> p and pi+; L = 1
        EvtResonance2 LambdacpKpi3(
            p0, p4Daug[2], p4Daug[0], 1.0, 0.0, _LambdaW, _LambdaM, 2, true,
            _LambdaR, _Lambda_cR );    // Lambda(1520) -> K- and p; L = 2

        std::vector<EvtComplex> ampNonRes =
            calcResAmpTerms( EvtLambdacPHH::NonReson, LambdacpKpi0, _NRNorm );
        std::vector<EvtComplex> ampKstar =
            calcResAmpTerms( EvtLambdacPHH::Kstar, LambdacpKpi1, _KstarNorm );
        std::vector<EvtComplex> ampDelta =
            calcResAmpTerms( EvtLambdacPHH::Delta, LambdacpKpi2, _DeltaNorm );
        std::vector<EvtComplex> ampLambda =
            calcResAmpTerms( EvtLambdacPHH::Lambda, LambdacpKpi3, _LambdaNorm );

        // Combine resonance amplitudes for a given spin configuration
        EvtComplex amp00 = ampNonRes[0] + ampKstar[0] + ampDelta[0] +
                           ampLambda[0];
        EvtComplex amp01 = ampNonRes[1] + ampKstar[1] + ampDelta[1] +
                           ampLambda[1];
        EvtComplex amp10 = ampNonRes[2] + ampKstar[2] + ampDelta[2] +
                           ampLambda[2];
        EvtComplex amp11 = ampNonRes[3] + ampKstar[3] + ampDelta[3] +
                           ampLambda[3];

        // Fit fraction numerator terms
        FitFracTop[0] += abs2( ampNonRes[0] ) + abs2( ampNonRes[1] ) +
                         abs2( ampNonRes[2] ) + abs2( ampNonRes[3] );
        FitFracTop[1] += abs2( ampKstar[0] ) + abs2( ampKstar[1] ) +
                         abs2( ampKstar[2] ) + abs2( ampKstar[3] );
        FitFracTop[2] += abs2( ampDelta[0] ) + abs2( ampDelta[1] ) +
                         abs2( ampDelta[2] ) + abs2( ampDelta[3] );
        FitFracTop[3] += abs2( ampLambda[0] ) + abs2( ampLambda[1] ) +
                         abs2( ampLambda[2] ) + abs2( ampLambda[3] );

        // Fit fraction common denominator
        FitFracDenom += abs2( amp00 ) + abs2( amp01 ) + abs2( amp10 ) +
                        abs2( amp11 );
    }

    EvtGenReport( EVTGEN_INFO, "EvtLambdacPHH" )
        << "FitFracs: NR = " << FitFracTop[0] / FitFracDenom
        << ", K* = " << FitFracTop[1] / FitFracDenom
        << ", Del = " << FitFracTop[2] / FitFracDenom
        << ", Lam = " << FitFracTop[3] / FitFracDenom << std::endl;
}

void EvtLambdacPHH::initProbMax()
{
    // Default value
    setProbMax( 10.0 );

    // Set probability using decay file parameter
    if ( getNArg() > 0 ) {
        setProbMax( getArg( 0 ) );
    }
}

void EvtLambdacPHH::decay( EvtParticle* p )
{
    // Daughter order: 1 = K-, 2 = pi+, 3 = p
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    // 4-momenta in the rest frame of the Lambda_c
    EvtVector4R p4_p( p->mass(), 0.0, 0.0, 0.0 );
    EvtVector4R moms1 = p->getDaug( _d1 )->getP4();
    EvtVector4R moms2 = p->getDaug( _d2 )->getP4();
    EvtVector4R moms3 = p->getDaug( _d3 )->getP4();

    // Lambda_c decay mode resonances. Spin L values from strong decay parity conservation:
    // parity(resonance) = parity(daug1)*parity(daug2)*(-1)^L
    EvtResonance2 LambdacpKpi0( p4_p, moms1, moms2, 1.0, 0.0, 0.0, 0.0, 0, true,
                                0.0, 0.0 );    // Non-resonant L = 0
    EvtResonance2 LambdacpKpi1( p4_p, moms1, moms2, 1.0, 0.0, _KstarW, _KstarM,
                                1, true, _KstarR,
                                _Lambda_cR );    // K*0 -> K- and pi+; L = 1
    EvtResonance2 LambdacpKpi2( p4_p, moms3, moms2, 1.0, 0.0, _DeltaW, _DeltaM,
                                1, true, _DeltaR,
                                _Lambda_cR );    // Delta++ -> p and pi+; L = 1
    EvtResonance2 LambdacpKpi3( p4_p, moms3, moms1, 1.0, 0.0, _LambdaW,
                                _LambdaM, 2, true, _LambdaR,
                                _Lambda_cR );    // Lambda(1520) -> K- and p; L = 2

    // Define the "beam" direction, used in Fig 1 of hep-ex/9912003v1
    EvtVector4R beam( 0.0, 0.0, 0.0, 1.0 );
    EvtParticle* parent = p->getParent();
    if ( parent ) {
        // If non prompt, the beam is along the direction of the mother
        EvtVector4R p4_Lambda_c_mother = parent->getP4Lab();
        p4_Lambda_c_mother.applyBoostTo( p->getP4Lab() );
        beam = p4_Lambda_c_mother;
    }

    _p4_Lambda_c = p->getP4Lab();
    _p4_Lambdac_Mag = _p4_Lambda_c.d3mag();

    // Define the unit vector denoting the "z" axis in Fig 1
    _zprime = -1.0 * _p4_Lambda_c.cross( beam );
    _zprime.applyBoostTo( _p4_Lambda_c, true );    // From lab frame to Lambda_c

    _zpMag = _zprime.d3mag();
    // Check if zprime magnitude is non-zero
    if ( _zpMag > 0.0 ) {
        // Normalise
        _zprime /= _zpMag;
    } else {
        // Set as the z direction
        _zprime.set( 0.0, 0.0, 0.0, 1.0 );
    }
    // Update normalised |z'|
    _zpMag = 1.0;

    // Get the amplitudes: non-resonant, K*, Delta and Lambda
    std::vector<EvtComplex> ampNonRes = calcResAmpTerms( EvtLambdacPHH::NonReson,
                                                         LambdacpKpi0, _NRNorm );
    std::vector<EvtComplex> ampKstar =
        calcResAmpTerms( EvtLambdacPHH::Kstar, LambdacpKpi1, _KstarNorm );
    std::vector<EvtComplex> ampDelta =
        calcResAmpTerms( EvtLambdacPHH::Delta, LambdacpKpi2, _DeltaNorm );
    std::vector<EvtComplex> ampLambda =
        calcResAmpTerms( EvtLambdacPHH::Lambda, LambdacpKpi3, _LambdaNorm );

    // Combine resonance amplitudes for a given spin configuration
    EvtComplex amp00 = ampNonRes[0] + ampKstar[0] + ampDelta[0] + ampLambda[0];
    EvtComplex amp01 = ampNonRes[1] + ampKstar[1] + ampDelta[1] + ampLambda[1];
    EvtComplex amp10 = ampNonRes[2] + ampKstar[2] + ampDelta[2] + ampLambda[2];
    EvtComplex amp11 = ampNonRes[3] + ampKstar[3] + ampDelta[3] + ampLambda[3];

    // Set the amplitude components
    vertex( 0, 0, amp00 );
    vertex( 0, 1, amp01 );
    vertex( 1, 0, amp10 );
    vertex( 1, 1, amp11 );
}

std::vector<EvtComplex> EvtLambdacPHH::calcResAmpTerms(
    EvtLambdacPHH::LcResLabel resIndex, const EvtResonance2& res, double norm ) const
{
    // Initialise the resonance and daughter theta and phi angles
    double thetaRes( 0.0 ), phiRes( 0.0 ), phiPrimeDaug( 0.0 ),
        thetaPrimeDaug( 0.0 );
    // Initialise beta rotation angle
    double beta_res( 0.0 );

    EvtVector4R res_atproton( 0.0, 0.0, 0.0, 0.0 ),
        Lc_atproton( 0.0, 0.0, 0.0, 0.0 );

    // Initialise Amplitude terms
    EvtComplex term1( 0.0 ), term2( 0.0 ), term3( 0.0 ), term4( 0.0 );
    // Normalised dynamical amplitude
    EvtComplex resAmp( norm, 0.0 );

    // Angles are not needed for the non-resonant amplitude
    if ( resIndex != EvtLambdacPHH::NonReson ) {
        resAmp = res.resAmpl() * norm;
        // Resonance and daughter 4 momenta
        EvtVector4R p4d1 = res.p4_d1();
        EvtVector4R p4d2 = res.p4_d2();
        EvtVector4R p4Res = p4d1 + p4d2;
        EvtVector4R p4_d3 = res.p4_p() - p4Res;

        double p4ResMag = p4Res.d3mag();

        // 4-momenta for theta' and phi' angles
        EvtVector4R yRes = -1.0 * p4_d3.cross( _zprime );

        EvtVector4R res_d1 = p4d1;
        res_d1.applyBoostTo( p4Res, true );
        double res_d1_Mag = res_d1.d3mag();

        EvtVector4R res_d3 = -1.0 * p4_d3;
        double res_d3_Mag = res_d3.d3mag();

        thetaPrimeDaug = getACos( res_d1.dot( res_d3 ), res_d1_Mag * res_d3_Mag );

        res_atproton = p4Res;
        res_atproton.applyBoostTo( p4d1, true );
        double res_atproton_mag = res_atproton.d3mag();

        Lc_atproton = res.p4_p();
        Lc_atproton.applyBoostTo( p4d1, true );
        double Lc_atproton_mag = Lc_atproton.d3mag();

        // Check that the momentum of the Lambda_c is not zero, as well as a valid zprime vector
        if ( _p4_Lambdac_Mag > 0.0 && _zpMag > 0.0 ) {
            thetaRes = getACos( -1.0 * p4Res.dot( _zprime ), p4ResMag );
            phiRes = getASin( -1.0 * p4Res.dot( _p4_Lambda_c ),
                              sin( thetaRes ) * _p4_Lambdac_Mag * p4ResMag );
            phiPrimeDaug = getASin( res_d1.dot( yRes ), sin( thetaPrimeDaug ) *
                                                            res_d1_Mag *
                                                            yRes.d3mag() );

        } else {
            // Use randomised angles with flat probability distributions
            thetaRes = EvtRandom::Flat( 0.0, EvtConst::pi );
            phiRes = EvtRandom::Flat( 0.0, EvtConst::twoPi );
            phiPrimeDaug = EvtRandom::Flat( 0.0, EvtConst::twoPi );
        }

        if ( res_atproton_mag > 0.0 && Lc_atproton_mag > 0.0 ) {
            // Extra rotation to go to the proton helicity frame for the two resonances Delta++ and Lambda.
            // No rotation is needed for K*. Use the momenta boosted to the proton restframe

            beta_res = getACos( res_atproton.dot( Lc_atproton ),
                                res_atproton_mag * Lc_atproton_mag );

        } else {
            beta_res = EvtRandom::Flat( 0.0, EvtConst::pi );
        }
    }

    // Find the spin-dependent amplitudes
    if ( resIndex == EvtLambdacPHH::NonReson ||
         resIndex == EvtLambdacPHH::Kstar ) {
        term1 = resAmp * DecayAmp3( resIndex, 1, 1, thetaRes, phiRes,
                                    thetaPrimeDaug, phiPrimeDaug );
        term2 = resAmp * DecayAmp3( resIndex, 1, -1, thetaRes, phiRes,
                                    thetaPrimeDaug, phiPrimeDaug );
        term3 = resAmp * DecayAmp3( resIndex, -1, 1, thetaRes, phiRes,
                                    thetaPrimeDaug, phiPrimeDaug );
        term4 = resAmp * DecayAmp3( resIndex, -1, -1, thetaRes, phiRes,
                                    thetaPrimeDaug, phiPrimeDaug );

    } else {
        double rotate_00 = EvtdFunction::d( 1, 1, 1, beta_res );
        double rotate_10 = EvtdFunction::d( 1, -1, 1, beta_res );
        double rotate_11 = EvtdFunction::d( 1, -1, -1, beta_res );
        double rotate_01 = EvtdFunction::d( 1, 1, -1, beta_res );

        // Delta and Lambda need to be rotated before summing over the proton helicity axis
        EvtComplex termA = resAmp * DecayAmp3( resIndex, 1, 1, thetaRes, phiRes,
                                               thetaPrimeDaug, phiPrimeDaug );
        EvtComplex termB = resAmp * DecayAmp3( resIndex, 1, -1, thetaRes, phiRes,
                                               thetaPrimeDaug, phiPrimeDaug );
        EvtComplex termC = resAmp * DecayAmp3( resIndex, -1, 1, thetaRes, phiRes,
                                               thetaPrimeDaug, phiPrimeDaug );
        EvtComplex termD = resAmp * DecayAmp3( resIndex, -1, -1, thetaRes, phiRes,
                                               thetaPrimeDaug, phiPrimeDaug );

        term1 = rotate_00 * termA + rotate_10 * termB;
        term2 = rotate_01 * termA + rotate_11 * termB;
        term3 = rotate_00 * termC + rotate_10 * termD;
        term4 = rotate_01 * termC + rotate_11 * termD;
    }

    // Return the spin amplitudes as a vector
    std::vector<EvtComplex> ampVect;
    ampVect.push_back( term1 );
    ampVect.push_back( term2 );
    ampVect.push_back( term3 );
    ampVect.push_back( term4 );

    return ampVect;
}

EvtComplex EvtLambdacPHH::DecayAmp3( EvtLambdacPHH::LcResLabel resonance, int m,
                                     int mprime, double theta_res, double phi_res,
                                     double theta_prime_daughter_res,
                                     double phi_prime_daughter_res ) const
{
    // Find the amplitudes given in Tables 3 to 6 in the paper.
    // Wigner d-functions use 2*spin, e.g. d(1/2, 1/2, 1/2) -> d(1, 1, 1)
    EvtComplex term1( 0.0, 0.0 ), term2( 0.0, 0.0 );

    if ( resonance == EvtLambdacPHH::NonReson ) {
        // Non-resonant: table 6
        if ( m == 1 && mprime == 1 ) {
            term1 = _Nplusplus *
                    EvtComplex( cos( _phiNplusplus ), sin( _phiNplusplus ) );

        } else if ( m == 1 && mprime == -1 ) {
            term1 = _Nplusminus *
                    EvtComplex( cos( _phiNplusminus ), sin( _phiNplusminus ) );

        } else if ( m == -1 && mprime == 1 ) {
            term1 = _Nminusplus *
                    EvtComplex( cos( _phiNminusplus ), sin( _phiNminusplus ) );

        } else if ( m == -1 && mprime == -1 ) {
            term1 = _Nminusminus * EvtComplex( cos( _phiNminusminus ),
                                               sin( _phiNminusminus ) );
        }

    } else if ( resonance == EvtLambdacPHH::Kstar ) {
        // K*0(1-) resonance: table 3
        if ( m == 1 && mprime == 1 ) {
            term1 = fampl3( _E1, _phiE1, 1, 1, 1, theta_res, 2, 2, 0,
                            theta_prime_daughter_res, phi_prime_daughter_res );
            term2 = fampl3( _E2, _phiE2, 1, 1, -1, theta_res, 2, 0, 0,
                            theta_prime_daughter_res, phi_res );

        } else if ( m == 1 && mprime == -1 ) {
            term1 = fampl3( _E3, _phiE3, 1, 1, 1, theta_res, 2, 0, 0,
                            theta_prime_daughter_res, 0.0 );
            term2 = fampl3( _E4, _phiE4, 1, 1, -1, theta_res, 2, -2, 0,
                            theta_prime_daughter_res,
                            phi_res - phi_prime_daughter_res );

        } else if ( m == -1 && mprime == 1 ) {
            term1 = fampl3( _E1, _phiE1, 1, -1, 1, theta_res, 2, 2, 0,
                            theta_prime_daughter_res,
                            -( phi_res - phi_prime_daughter_res ) );
            term2 = fampl3( _E2, _phiE2, 1, -1, -1, theta_res, 2, 0, 0,
                            theta_prime_daughter_res, 0.0 );

        } else if ( m == -1 && mprime == -1 ) {
            term1 = fampl3( _E3, _phiE3, 1, -1, 1, theta_res, 2, 0, 0,
                            theta_prime_daughter_res, -phi_res );
            term2 = fampl3( _E4, _phiE4, 1, -1, -1, theta_res, 2, -2, 0,
                            theta_prime_daughter_res, -phi_prime_daughter_res );
        }

    } else if ( resonance == EvtLambdacPHH::Delta ) {
        // Delta++(3/2+) resonance: table 4
        if ( m == 1 && mprime == 1 ) {
            term1 = fampl3( _F1, _phiF1, 1, 1, 1, theta_res, 3, 1, 1,
                            theta_prime_daughter_res, 0.0 );
            term2 = fampl3( _F2, _phiF2, 1, 1, -1, theta_res, 3, -1, 1,
                            theta_prime_daughter_res,
                            phi_res - phi_prime_daughter_res );

        } else if ( m == 1 && mprime == -1 ) {
            term1 = fampl3( _F1, _phiF1, 1, 1, 1, theta_res, 3, 1, -1,
                            theta_prime_daughter_res, phi_prime_daughter_res );
            term2 = fampl3( _F2, _phiF2, 1, 1, -1, theta_res, 3, -1, -1,
                            theta_prime_daughter_res, phi_res );

        } else if ( m == -1 && mprime == 1 ) {
            term1 = fampl3( _F1, _phiF1, 1, -1, 1, theta_res, 3, 1, 1,
                            theta_prime_daughter_res, -phi_res );
            term2 = fampl3( _F2, _phiF2, 1, -1, -1, theta_res, 3, -1, 1,
                            theta_prime_daughter_res, -phi_prime_daughter_res );

        } else if ( m == -1 && mprime == -1 ) {
            term1 = fampl3( _F1, _phiF1, 1, -1, 1, theta_res, 3, 1, -1,
                            theta_prime_daughter_res,
                            -( phi_res - phi_prime_daughter_res ) );
            term2 = fampl3( _F2, _phiF2, 1, -1, -1, theta_res, 3, -1, -1,
                            theta_prime_daughter_res, 0.0 );
        }

    } else if ( resonance == EvtLambdacPHH::Lambda ) {
        // Lambda(1520)(3/2-) resonance: table 5
        if ( m == 1 && mprime == 1 ) {
            term1 = fampl3( _H1, _phiH1, 1, 1, 1, theta_res, 3, 1, 1,
                            theta_prime_daughter_res, 0.0 );
            term2 = fampl3( _H2, _phiH2, 1, 1, -1, theta_res, 3, -1, 1,
                            theta_prime_daughter_res,
                            phi_res - phi_prime_daughter_res );

        } else if ( m == 1 && mprime == -1 ) {
            term1 = -1.0 * fampl3( _H1, _phiH1, 1, 1, 1, theta_res, 3, 1, -1,
                                   theta_prime_daughter_res,
                                   phi_prime_daughter_res );
            term2 = -1.0 * fampl3( _H2, _phiH2, 1, 1, -1, theta_res, 3, -1, -1,
                                   theta_prime_daughter_res, phi_res );

        } else if ( m == -1 && mprime == 1 ) {
            term1 = fampl3( _H1, _phiH1, 1, -1, 1, theta_res, 3, 1, 1,
                            theta_prime_daughter_res, -phi_res );
            term2 = fampl3( _H2, _phiH2, 1, -1, -1, theta_res, 3, -1, 1,
                            theta_prime_daughter_res, -phi_prime_daughter_res );

        } else if ( m == -1 && mprime == -1 ) {
            term1 = -1.0 * fampl3( _H1, _phiH1, 1, -1, 1, theta_res, 3, 1, -1,
                                   theta_prime_daughter_res,
                                   -( phi_res - phi_prime_daughter_res ) );
            term2 = -1.0 * fampl3( _H2, _phiH2, 1, -1, -1, theta_res, 3, -1, -1,
                                   theta_prime_daughter_res, 0.0 );
        }
    }

    EvtComplex Amplitude = term1 + term2;
    return Amplitude;
}

EvtComplex EvtLambdacPHH::fampl3( double amplitude_res, double phi_res,
                                  int spinMother, int m_spinMother,
                                  int m_prime_spinMother, double theta_res,
                                  float spin_res, float m_spin_res,
                                  float m_prime_spin_res,
                                  double theta_daughter_res,
                                  double phi_prime_daughter_res ) const
{
    double dTerm1 = EvtdFunction::d( spinMother, m_spinMother,
                                     m_prime_spinMother, theta_res );
    double dTerm2 = EvtdFunction::d( spin_res, m_spin_res, m_prime_spin_res,
                                     theta_daughter_res );

    EvtComplex amp_phase1 = EvtComplex( cos( phi_res ), sin( phi_res ) );
    EvtComplex amp_phase2 = EvtComplex( cos( phi_prime_daughter_res ),
                                        sin( phi_prime_daughter_res ) );

    EvtComplex partial_amp = amplitude_res * amp_phase1 * dTerm1 * amp_phase2 *
                             dTerm2;

    return partial_amp;
}

double EvtLambdacPHH::getACos( double num, double denom ) const
{
    // Find inverse cosine, checking ratio is within +- 1
    double angle( 0.0 ), ratio( 0.0 );
    if ( fabs( denom ) > 0.0 ) {
        ratio = num / denom;
    }

    if ( fabs( ratio ) <= 1.0 ) {
        angle = acos( ratio );
    }

    return angle;
}

double EvtLambdacPHH::getASin( double num, double denom ) const
{
    // Find inverse sine, checking ratio is within +- 1
    double angle( 0.0 ), ratio( 0.0 );
    if ( fabs( denom ) > 0.0 ) {
        ratio = num / denom;
    }

    if ( fabs( ratio ) <= 1.0 ) {
        angle = asin( ratio );
    }

    return angle;
}
