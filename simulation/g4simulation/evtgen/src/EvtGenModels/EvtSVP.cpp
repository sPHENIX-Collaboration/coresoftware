
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

#include "EvtGenModels/EvtSVP.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <cmath>

std::string EvtSVP::getName()
{
    return "SVP";
}

EvtDecayBase* EvtSVP::clone()
{
    return new EvtSVP;
}

void EvtSVP::decay_2body( EvtParticle* root )
{
    root->initializePhaseSpace( getNDaug(), getDaugs() );
    // Photon is the first particle and psi is the second
    // to ensure decay file backwards compatibility
    EvtParticle* photon = root->getDaug( 0 );
    EvtParticle* psi = root->getDaug( 1 );

    EvtVector4R p = psi->getP4(),    // psi momentum
        k = photon->getP4();         // Photon momentum

    bool validAmp( true );

    double kp = k * p;
    if ( fabs( kp ) < 1e-10 ) {
        validAmp = false;
    }

    for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
        EvtVector4C epsPsi = psi->epsParent( iPsi ).conj();
        for ( int iGamma = 0; iGamma < 2; iGamma++ ) {
            EvtVector4C epsGamma = photon->epsParentPhoton( iGamma ).conj();
            EvtComplex amp( 0.0, 0.0 );
            if ( validAmp ) {
                amp = ( epsPsi * epsGamma ) -
                      ( epsPsi * k ) * ( epsGamma * p ) / kp;
            }
            vertex( iGamma, iPsi, amp );
        }
    }
}

void EvtSVP::decay_3body( EvtParticle* root )
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

    // Extra checks to make sure we are not dividing by zero
    double kp = k * p;
    if ( fabs( kp ) < 1e-10 ) {
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
    for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
        EvtVector4C epsPsi = psi->epsParent( iPsi ).conj();

        for ( int iMplus = 0; iMplus < 2; iMplus++ ) {
            EvtDiracSpinor spMplus = mup->spParent( iMplus );

            for ( int iMminus = 0; iMminus < 2; iMminus++ ) {
                EvtDiracSpinor spMminus = mum->spParent( iMminus );
                EvtVector4C epsGamma = EvtLeptonVCurrent( spMplus, spMminus );
                EvtComplex amp( 0.0, 0.0 );
                if ( validAmp ) {
                    amp = ( epsPsi * epsGamma ) -
                          ( epsPsi * k ) * ( epsGamma * p ) / kp;
                }
                amp *= factor;

                // Set the amplitude matrix element using the vertex function
                vertex( iPsi, iMplus, iMminus, amp );
            }
        }
    }
}

void EvtSVP::decay( EvtParticle* root )
{
    if ( getNDaug() == 2 ) {
        decay_2body( root );
    } else if ( getNDaug() == 3 ) {
        decay_3body( root );
    }
}

void EvtSVP::init()
{
    checkSpinParent( EvtSpinType::SCALAR );

    if ( getNDaug() == 2 ) {    // chi -> gamma psi radiative mode
        checkNArg( 0 );
        checkNDaug( 2 );
        checkSpinDaughter( 0, EvtSpinType::PHOTON );
        checkSpinDaughter( 1, EvtSpinType::VECTOR );

    } else if ( getNDaug() == 3 ) {    // chi -> psi lepton lepton
        checkSpinParent( EvtSpinType::SCALAR );
        checkSpinDaughter( 0, EvtSpinType::VECTOR );
        checkSpinDaughter( 1, EvtSpinType::DIRAC );
        checkSpinDaughter( 2, EvtSpinType::DIRAC );
        checkNArg( 1 );
        delta = getArg( 0 );
    }
}

void EvtSVP::initProbMax()
{
    if ( getNDaug() == 2 ) {
        setProbMax( 2.2 );

    } else if ( getNDaug() == 3 ) {
        const EvtId daugId = getDaug( 1 );

        if ( daugId == EvtPDL::getId( "mu+" ) ||
             daugId == EvtPDL::getId( "mu-" ) ) {
            setProbMax( 130.0 );
        } else if ( daugId == EvtPDL::getId( "e+" ) ||
                    daugId == EvtPDL::getId( "e-" ) ) {
            setProbMax( 4100.0 );
        }
    }
}
