
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

#include "EvtGenModels/EvtEtaLLPiPi.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <cmath>

// eta' -> mu+ mu- pi+ pi- or e+ e- pi+ pi-
// From Zhang Zhen-Yu et al, Chinese Phys. C 36, p926, 2012

void EvtEtaLLPiPi::init()
{
    // Check for 0 or 1 (maxProb) arguments
    checkNArg( 0, 1 );

    // Check particle types
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::SCALAR );
    checkSpinDaughter( 3, EvtSpinType::SCALAR );

    // Mass and width of rho0 from particle properties file
    m_rhoMass = EvtPDL::getMeanMass( EvtPDL::getId( "rho0" ) );
    m_rhoMassSq = m_rhoMass * m_rhoMass;
    m_rhoGamma = EvtPDL::getWidth( EvtPDL::getId( "rho0" ) );

    // Mixing parameter squared, using Eq 6
    const double denom = 8.0 * pow( EvtConst::pi * m_fPi, 2 );
    const double factor = m_eSq / ( denom * denom * 3.0 );
    const double fTerm8 = sin( m_thetaMix ) / m_f8;
    const double fTerm0 = 2.0 * sqrt( 2.0 ) * cos( m_thetaMix ) / m_f0;
    m_mixSq = factor * pow( fTerm8 + fTerm0, 2 );
}

void EvtEtaLLPiPi::initProbMax()
{
    if ( getNArg() == 1 ) {
        setProbMax( getArg( 0 ) );

    } else {
        int lepId = getDaug( 0 ).getId();
        if ( lepId == EvtPDL::getId( "e-" ).getId() ||
             lepId == EvtPDL::getId( "e+" ).getId() ) {
            setProbMax( 27500.0 );

        } else if ( lepId == EvtPDL::getId( "mu-" ).getId() ||
                    lepId == EvtPDL::getId( "mu+" ).getId() ) {
            setProbMax( 20.0 );
        }
    }
}

std::string EvtEtaLLPiPi::getName()
{
    return "ETA_LLPIPI";
}

EvtDecayBase* EvtEtaLLPiPi::clone()
{
    return new EvtEtaLLPiPi();
}

void EvtEtaLLPiPi::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    const double mLep = p->getDaug( 0 )->mass();
    const double mPi = p->getDaug( 2 )->mass();

    updateMassPars( mLep, mPi );

    const double prob = ampSquared( p );
    setProb( prob );
}

void EvtEtaLLPiPi::updateMassPars( double mLep, double mPi )
{
    // Update mass parameters used in various functions
    m_lepMass = mLep;
    m_lepMassSq = mLep * mLep;
    m_4LepMassSq = 4.0 * m_lepMassSq;

    m_piMass = mPi;
    m_piMassSq = mPi * mPi;
    m_4PiMassSq = 4.0 * m_piMassSq;
}

double EvtEtaLLPiPi::rhoWidth( double s, double m ) const
{
    // Define width of rho using modified vector meson dynamics, Eq 8
    double gamma( 0.0 );

    if ( s >= m_4PiMassSq ) {
        const double mSq = m * m;
        const double num = 1.0 - ( 4.0 * mSq / s );
        const double denom = 1.0 - ( 4.0 * mSq / m_rhoMassSq );
        const double ratio = denom > 0.0 ? num / denom : 0.0;
        gamma = m_rhoGamma * ( s / m_rhoMassSq ) * pow( ratio, 1.5 );
    }

    return gamma;
}

double EvtEtaLLPiPi::F0( double sLL, double sPiPi ) const
{
    // Equation 7
    double ampSq( 0.0 );

    const double rhoWidthL = rhoWidth( sLL, m_lepMass );
    const double rhoWidthPi = rhoWidth( sPiPi, m_piMass );

    const double mSqDiffL = m_rhoMassSq - sLL;
    const double mRhoWidthL = m_rhoMass * rhoWidthL;

    const double mSqDiffPi = m_rhoMassSq - sPiPi;
    const double mRhoWidthPi = m_rhoMass * rhoWidthPi;

    const double denomLL = mSqDiffL * mSqDiffL + mRhoWidthL * mRhoWidthL;
    const double denomPiPi = mSqDiffPi * mSqDiffPi + mRhoWidthPi * mRhoWidthPi;

    if ( denomLL > 0.0 && denomPiPi > 0.0 ) {
        const double mRho4 = m_rhoMassSq * m_rhoMassSq;
        const double denomProd = denomLL * denomPiPi;

        double realAmp = m_par1 + m_parLL * ( m_rhoMassSq * mSqDiffL / denomLL );
        realAmp += m_parPiPi * mRho4 *
                   ( ( mSqDiffPi * mSqDiffL ) - mRhoWidthL * mRhoWidthPi ) /
                   denomProd;

        double imagAmp = m_parLL * ( m_rhoMassSq * mRhoWidthL / denomLL );
        imagAmp += m_parPiPi * mRho4 *
                   ( mRhoWidthPi * mSqDiffL + mRhoWidthL * mSqDiffPi ) /
                   denomProd;

        ampSq = realAmp * realAmp + imagAmp * imagAmp;
    }

    return ampSq;
}

double EvtEtaLLPiPi::lambda( double a, double b, double c ) const
{
    const double sumSq = a * a + b * b + c * c;
    const double prod = a * b + b * c + c * a;
    const double L = sumSq - 2.0 * prod;
    return L;
}

double EvtEtaLLPiPi::ampSquared( EvtParticle* p ) const
{
    // Equation 3
    const double zeroProb( 0.0 );

    // Mass of eta' meson
    const double mEta = p->mass();

    const EvtVector4R pL1 = p->getDaug( 0 )->getP4();
    const EvtVector4R pL2 = p->getDaug( 1 )->getP4();
    const EvtVector4R pPi1 = p->getDaug( 2 )->getP4();
    const EvtVector4R pPi2 = p->getDaug( 3 )->getP4();

    const EvtVector4R pLL = pL1 + pL2;
    const double sLL = pLL.mass2();
    const EvtVector4R pPiPi = pPi1 + pPi2;
    const double sPiPi = pPiPi.mass2();

    if ( sLL < 1e-4 || sPiPi < m_4PiMassSq || sLL < m_4LepMassSq ) {
        // To avoid negative square roots etc
        return zeroProb;
    }

    // Angles theta_p, theta_k and phi defined in Fig 1
    const EvtVector4R pSum = pLL + pPiPi;
    // Helicity angle of first lepton
    const double cosThp = -EvtDecayAngle( pSum, pLL, pL1 );
    const double sinThp = sqrt( 1.0 - cosThp * cosThp );
    // Helicity angle of first pion
    const double cosThk = -EvtDecayAngle( pSum, pPiPi, pPi2 );
    const double sinThk = sqrt( 1.0 - cosThk * cosThk );
    // Angle between the dilepton and dipion decay planes
    const double phi = EvtDecayAngleChi( pSum, pL1, pL2, pPi1, pPi2 );
    const double sinPhi = sin( phi );

    const double betaLL = sqrt( 1.0 - ( m_4LepMassSq / sLL ) );
    const double betaPiPi = sqrt( 1.0 - ( m_4PiMassSq / sPiPi ) );

    const double betaProd = ( 1.0 - pow( betaLL * sinThp * sinPhi, 2 ) ) *
                            sPiPi * pow( betaPiPi * sinThk, 2 );
    const double L = lambda( mEta * mEta, sLL, sPiPi );
    const double ampSq = m_eSq * F0( sLL, sPiPi ) * m_mixSq * L * betaProd /
                         ( 8.0 * sLL );

    return ampSq;
}
