
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

#include "EvtGenModels/EvtVVP.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <cmath>

std::string EvtVVP::getName()
{
    return "VVP";
}

EvtDecayBase* EvtVVP::clone()
{
    return new EvtVVP;
}

void EvtVVP::init()
{
    checkSpinParent( EvtSpinType::VECTOR );

    if ( getNDaug() == 2 ) {    // chi -> gamma psi radiative mode
        // This model needs 0 parameters, but previously was defined as requiring 8!
        // Check for 0 or 8 parameters in the decay file for backwards compatibility
        checkNArg( 0, 8 );
        checkNDaug( 2 );
        checkSpinDaughter( 0, EvtSpinType::VECTOR );
        checkSpinDaughter( 1, EvtSpinType::PHOTON );

    } else if ( getNDaug() == 3 ) {    // chi -> psi lepton lepton
        checkSpinDaughter( 0, EvtSpinType::VECTOR );
        checkSpinDaughter( 1, EvtSpinType::DIRAC );
        checkSpinDaughter( 2, EvtSpinType::DIRAC );
        checkNArg( 1 );
        delta = getArg( 0 );
    }
}

void EvtVVP::initProbMax()
{
    if ( getNDaug() == 2 ) {
        setProbMax( 2.0 );

    } else if ( getNDaug() == 3 ) {
        const EvtId daugId = getDaug( 1 );

        if ( daugId == EvtPDL::getId( "mu+" ) ||
             daugId == EvtPDL::getId( "mu-" ) ) {
            setProbMax( 15.0 );
        } else if ( daugId == EvtPDL::getId( "e+" ) ||
                    daugId == EvtPDL::getId( "e-" ) ) {
            setProbMax( 600.0 );
        }
    }
}

void EvtVVP::decay( EvtParticle* root )
{
    if ( getNDaug() == 2 ) {
        decay_2body( root );
    } else if ( getNDaug() == 3 ) {
        decay_3body( root );
    }
}

void EvtVVP::decay_2body( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    // Vector is first particle, photon is the second
    EvtParticle *v, *ph;
    v = p->getDaug( 0 );
    ph = p->getDaug( 1 );
    EvtVector3C epsp[3];
    EvtVector3C epsv[3];
    EvtVector3C epsph[2];
    epsp[0] = p->eps( 0 ).vec();
    epsp[1] = p->eps( 1 ).vec();
    epsp[2] = p->eps( 2 ).vec();

    epsv[0] = v->eps( 0 ).vec().conj();
    epsv[1] = v->eps( 1 ).vec().conj();
    epsv[2] = v->eps( 2 ).vec().conj();

    epsph[0] = ph->epsParentPhoton( 0 ).vec().conj();
    epsph[1] = ph->epsParentPhoton( 1 ).vec().conj();

    int i, j, k;
    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < 3; j++ ) {
            for ( k = 0; k < 2; k++ ) {
                vertex( i, j, k, epsp[i].cross( epsv[j] ) * epsph[k] );
            }
        }
    }
}

void EvtVVP::decay_3body( EvtParticle* root )
{
    root->initializePhaseSpace( getNDaug(), getDaugs() );
    EvtParticle* psi = root->getDaug( 0 );
    EvtParticle* mup = root->getDaug( 1 );
    EvtParticle* mum = root->getDaug( 2 );

    EvtVector4R k1 = mup->getP4(),    // mu+ momentum
        k2 = mum->getP4(),            // mu- momentum
        k = k1 + k2;                  // photon momentum

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
    double dSq = delta * delta;
    double dSqDenom = dSq - kSq;
    if ( fabs( dSqDenom ) < 1e-10 ) {
        validAmp = false;
    }

    double factor( 1.0 );
    if ( validAmp ) {
        factor = dSq / ( dSqDenom * kSq );
    }

    int iPols[4] = {0, 0, 0, 0};

    // Calculate the amplitude terms, looping over the chi, psi and lepton states
    for ( int iChi = 0; iChi < 3; iChi++ ) {
        iPols[0] = iChi;
        EvtVector4C epsChi = root->epsParent( iChi );

        for ( int iPsi = 0; iPsi < 3; iPsi++ ) {
            iPols[1] = iPsi;
            EvtVector4C epsPsi = psi->epsParent( iPsi ).conj();

            for ( int iMplus = 0; iMplus < 2; iMplus++ ) {
                iPols[2] = iMplus;
                EvtDiracSpinor spMplus = mup->spParent( iMplus );

                for ( int iMminus = 0; iMminus < 2; iMminus++ ) {
                    iPols[3] = iMminus;
                    EvtDiracSpinor spMminus = mum->spParent( iMminus );
                    EvtVector4C epsGamma =
                        EvtLeptonVCurrent( spMplus, spMminus ).conj();

                    // Based on Baranov PRD 85,014034 (2012), Eq 10
                    // amp = e_{mu nu alpha beta} epsChi^mu epsPsi^nu epsGamma^alpha k^beta/k^2
                    EvtComplex amp( 0.0, 0.0 );
                    if ( validAmp ) {
                        amp = k * dual( EvtGenFunctions::directProd( epsChi,
                                                                     epsPsi ) )
                                      .cont1( epsGamma );
                    }
                    amp *= factor;

                    // Set the amplitude matrix element using the vertex function
                    vertex( iPols, amp );
                }
            }
        }
    }
}
