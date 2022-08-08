
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

#include "EvtGenModels/EvtTVP.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <cmath>

std::string EvtTVP::getName()
{
    return "TVP";
}

EvtDecayBase* EvtTVP::clone()
{
    return new EvtTVP;
}

void EvtTVP::decay( EvtParticle* root )
{
    if ( getNDaug() == 2 ) {
        decay_2body( root );
    } else if ( getNDaug() == 3 ) {
        decay_3body( root );
    }
}

void EvtTVP::init()
{
    checkSpinParent( EvtSpinType::TENSOR );

    if ( getNDaug() == 2 ) {    // chi -> gamma psi radiative mode
        checkNArg( 0 );
        checkSpinDaughter( 0, EvtSpinType::PHOTON );
        checkSpinDaughter( 1, EvtSpinType::VECTOR );
    } else if ( getNDaug() == 3 ) {    // chi -> psi lepton lepton
        checkNDaug( 3 );
        checkSpinDaughter( 0, EvtSpinType::VECTOR );
        checkSpinDaughter( 1, EvtSpinType::DIRAC );
        checkSpinDaughter( 2, EvtSpinType::DIRAC );
        checkNArg( 1 );
        delta = getArg( 0 );
    }
}

void EvtTVP::initProbMax()
{
    if ( getNDaug() == 2 ) {
        const EvtId parId = getParentId();
        if ( parId == EvtPDL::getId( "chi_b2" ) ) {
            setProbMax( 15.0 );
        } else {
            setProbMax( 2.0 );
        }

    } else if ( getNDaug() == 3 ) {
        double dSq = delta * delta;
        double denom = dSq - 0.2;
        double ratio( 1.0 );
        if ( fabs( denom ) > 1e-10 ) {
            ratio = dSq / denom;
        }
        double ffCor = ratio * ratio;

        const EvtId daugId = getDaug( 1 );
        const EvtId parId = getParentId();

        if ( daugId == EvtPDL::getId( "mu+" ) ||
             daugId == EvtPDL::getId( "mu-" ) ) {
            if ( parId == EvtPDL::getId( "chi_c2" ) ) {
                setProbMax( ffCor * 85.0 );    // tested on 1e6 events
            } else if ( parId == EvtPDL::getId( "chi_b2" ) ) {
                setProbMax( ffCor * 750.0 );    // tested on 1e6 events
            }

        } else if ( daugId == EvtPDL::getId( "e+" ) ||
                    daugId == EvtPDL::getId( "e-" ) ) {
            if ( parId == EvtPDL::getId( "chi_c2" ) ) {
                setProbMax( ffCor * 3.5e3 );    // tested on 1e5 events
            } else if ( parId == EvtPDL::getId( "chi_b2" ) ) {
                setProbMax( ffCor * 2.6e4 );
            }
        }
    }
}

void EvtTVP::decay_2body( EvtParticle* root )
{
    root->initializePhaseSpace( getNDaug(), getDaugs() );

    // Photon is the first particle and psi is the second
    // to ensure decay file backwards compatibility
    EvtParticle* photon = root->getDaug( 0 );
    EvtParticle* psi = root->getDaug( 1 );

    EvtVector4R p = psi->getP4(),    // psi momentum
        k = photon->getP4();         // Photon momentum

    for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
        EvtVector4C epsPsi = psi->epsParent( iPsi ).conj();

        for ( int iGamma = 0; iGamma < 2; iGamma++ ) {
            EvtVector4C epsGamma = photon->epsParentPhoton( iGamma ).conj();

            for ( int iChi = 0; iChi < 5; iChi++ ) {
                EvtTensor4C epsChi = root->epsTensor( iChi );

                // Baranov PRD 85,014034 (2012), Eq 11
                // amp = p^mu epsPsi^a epsChi_{a b} [k_mu epsGamma_b  - k_b epsGamma_mu]
                EvtVector4C eee = epsChi.cont1( epsPsi );
                EvtVector4C vvv = ( p * k ) * eee - ( k * eee ) * p;
                EvtComplex amp = vvv * epsGamma;
                vertex( iChi, iGamma, iPsi, amp );
            }
        }
    }
}

void EvtTVP::decay_3body( EvtParticle* root )
{
    root->initializePhaseSpace( getNDaug(), getDaugs() );
    EvtParticle* psi = root->getDaug( 0 );
    EvtParticle* mup = root->getDaug( 1 );
    EvtParticle* mum = root->getDaug( 2 );

    EvtVector4R p = psi->getP4(),    // psi momentum
        k1 = mup->getP4(),           // mu+ momentum
        k2 = mum->getP4(),           // mu- momentum
        k = k1 + k2;                 // photon momentum

    double kSq = k * k;

    // The decay amplitude needs four-vector products. Make sure we have
    // valid values for these, otherwise set the amplitude to zero.
    // We need to set _amp2 (EvtDecayAmp) via the vertex() function call
    // even when the amplitude is zero, otherwise the amplitude from the
    // previous accepted event will be used, potentially leading to biases

    // Selection on k^2 to avoid inefficient generation for the electron modes
    bool validAmp( true );
    if ( kSq < 1e-3 ) {
        validAmp = false;
    }

    double dSq = delta * delta;
    double dSqDenom = dSq - kSq;
    if ( fabs( dSqDenom ) < 1e-10 ) {
        validAmp = false;
    }

    double factor( 1.0 );
    if ( validAmp ) {
        factor = dSq / ( dSqDenom * kSq );
    }

    // Calculate the amplitude terms, looping over the psi and lepton states
    int iPols[4] = {0, 0, 0, 0};

    for ( int iChi = 0; iChi < 5; iChi++ ) {
        iPols[0] = iChi;
        EvtTensor4C epsChi = root->epsTensor( iChi );

        for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
            iPols[1] = iPsi;
            EvtVector4C epsPsi = psi->epsParent( iPsi ).conj();

            for ( int iMplus = 0; iMplus < 2; iMplus++ ) {
                iPols[2] = iMplus;
                EvtDiracSpinor spMplus = mup->spParent( iMplus );

                for ( int iMminus = 0; iMminus < 2; iMminus++ ) {
                    iPols[3] = iMminus;
                    EvtDiracSpinor spMminus = mum->spParent( iMminus );
                    EvtVector4C epsGamma = EvtLeptonVCurrent( spMplus, spMminus );

                    // Based on Baranov PRD 85,014034 (2012), Eq 11
                    // amp = p^mu epsPsi^a epsChi_{a b} [k_mu epsGamma_b  - k_b epsGamma_mu]/k^2
                    EvtVector4C eee = epsChi.cont1( epsPsi );
                    EvtVector4C vvv = ( p * k ) * eee - ( k * eee ) * p;
                    EvtComplex amp( 0.0, 0.0 );
                    if ( validAmp ) {
                        amp = vvv * epsGamma;
                    }
                    amp *= factor;

                    // Set the amplitude matrix element using the vertex function
                    vertex( iPols, amp );
                }
            }
        }
    }
}
