
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

#include "EvtGenModels/EvtSLDiBaryonAmp.hh"

#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRaritaSchwinger.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"

EvtSLDiBaryonAmp::EvtSLDiBaryonAmp( const EvtBToDiBaryonlnupQCDFF& formFactors ) :
    ffModel_( formFactors )
{
}

void EvtSLDiBaryonAmp::CalcAmp( EvtParticle* parent, EvtAmp& amp ) const
{
    static EvtId EM = EvtPDL::getId( "e-" );
    static EvtId MUM = EvtPDL::getId( "mu-" );
    static EvtId TAUM = EvtPDL::getId( "tau-" );
    static EvtId EP = EvtPDL::getId( "e+" );
    static EvtId MUP = EvtPDL::getId( "mu+" );
    static EvtId TAUP = EvtPDL::getId( "tau+" );

    // The amplitude assumes B- -> p+ p- l- nubar ordering
    // i.e. the B- decay is the "particle" mode

    // B charge (x3) to check for antiparticle mode and baryon daughter ordering
    EvtId BId = parent->getId();
    int qB3 = EvtPDL::chg3( BId );

    bool particleMode( true );
    // Check if we have B+ instead (antiparticle mode)
    if ( qB3 > 0 ) {
        particleMode = false;
    }

    // The baryon, charged lepton and neutrino daughters

    // Make sure the first baryon has a charge opposite to the B, since the
    // amplitude expressions assume this order
    EvtParticle* baryon1 = parent->getDaug( 0 );
    EvtParticle* baryon2 = parent->getDaug( 1 );

    // Check if we need to reverse the baryon ordering
    if ( EvtPDL::chg3( baryon1->getId() ) == qB3 ) {
        baryon1 = parent->getDaug( 1 );
        baryon2 = parent->getDaug( 0 );
    }

    EvtParticle* lepton = parent->getDaug( 2 );
    EvtParticle* neutrino = parent->getDaug( 3 );

    // 4-momenta in B rest frame
    EvtVector4R p0( parent->mass(), 0.0, 0.0, 0.0 );
    EvtVector4R p1 = baryon1->getP4();
    EvtVector4R p2 = baryon2->getP4();

    EvtVector4R pSum = p1 + p2;
    EvtVector4R p = p0 - pSum;
    EvtVector4R pDiff = p2 - p1;

    // Particle id's: retrieve 1st baryon again in case order has changed
    EvtId Id1 = baryon1->getId();
    EvtId Id2 = baryon2->getId();
    EvtId l_num = lepton->getId();

    EvtSpinType::spintype type1 = EvtPDL::getSpinType( Id1 );
    EvtSpinType::spintype type2 = EvtPDL::getSpinType( Id2 );

    // Parity of B+- = -1. Check if the parity of the dibaryon state is the same.
    // If so, set the sameParity integer to 1. Otherwise set it to -1,
    // i.e. the dibaryon system has opposite parity to the B meson
    int J1 = EvtSpinType::getSpin2( type1 );
    int J2 = EvtSpinType::getSpin2( type2 );
    int sameParity = this->checkDibaryonParity( Id1, Id2, J1, J2 );

    // Number of chiral components of the baryon spinors
    int N1 = EvtSpinType::getSpinStates( type1 );
    int N2 = EvtSpinType::getSpinStates( type2 );

    // Invariant mass of the two baryon particle system
    double m_dibaryon = sqrt( pSum.mass2() );

    // Complex number i
    EvtComplex I( 0, 1 );

    // Lepton currents, same for all baryon options
    EvtVector4C l1, l2;

    if ( l_num == EM || l_num == MUM || l_num == TAUM ) {
        // B-
        l1 = EvtLeptonVACurrent( lepton->spParent( 0 ),
                                 neutrino->spParentNeutrino() );

        l2 = EvtLeptonVACurrent( lepton->spParent( 1 ),
                                 neutrino->spParentNeutrino() );

    } else if ( l_num == EP || l_num == MUP || l_num == TAUP ) {
        // B+
        l1 = EvtLeptonVACurrent( neutrino->spParentNeutrino(),
                                 lepton->spParent( 0 ) );

        l2 = EvtLeptonVACurrent( neutrino->spParentNeutrino(),
                                 lepton->spParent( 1 ) );

    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtSLDiBaryonAmp" )
            << "Wrong lepton number" << std::endl;
    }

    // Parity multiplication factors for the antiparticle mode hadronic currents
    double sign1 = ( particleMode == true ) ? 1.0 : 1.0 * sameParity;
    double sign2 = ( particleMode == true ) ? 1.0 : 1.0 * sameParity;
    double sign3 = ( particleMode == true ) ? 1.0 : -1.0 * sameParity;
    double sign4 = ( particleMode == true ) ? 1.0 : -1.0 * sameParity;
    double sign5 = ( particleMode == true ) ? 1.0 : -1.0 * sameParity;
    double sign6 = ( particleMode == true ) ? 1.0 : 1.0 * sameParity;

    // Define form factor coeff variables
    double f1( 0.0 ), f2( 0.0 ), f3( 0.0 ), f4( 0.0 ), f5( 0.0 );
    double g1( 0.0 ), g2( 0.0 ), g3( 0.0 ), g4( 0.0 ), g5( 0.0 );

    // Handle case of two Dirac-type daughters, e.g. p pbar, p N(1440)
    if ( type1 == EvtSpinType::DIRAC && type2 == EvtSpinType::DIRAC ) {
        // Form factor parameters
        EvtBToDiBaryonlnupQCDFF::FormFactors FF;
        ffModel_.getDiracFF( parent, m_dibaryon, FF );

        if ( sameParity == 1 ) {
            f1 = FF.F1;
            f2 = FF.F2;
            f3 = FF.F3;
            f4 = FF.F4;
            f5 = FF.F5;
            g1 = FF.G1;
            g2 = FF.G2;
            g3 = FF.G3;
            g4 = FF.G4;
            g5 = FF.G5;
        } else {
            // Swap coeffs: f_i <--> g_i
            f1 = FF.G1;
            f2 = FF.G2;
            f3 = FF.G3;
            f4 = FF.G4;
            f5 = FF.G5;
            g1 = FF.F1;
            g2 = FF.F2;
            g3 = FF.F3;
            g4 = FF.F4;
            g5 = FF.F5;
        }

        EvtVector4R gMtmTerms = g3 * p + g4 * pSum + g5 * pDiff;
        EvtVector4R fMtmTerms = f3 * p + f4 * pSum + f5 * pDiff;

        // First baryon
        for ( int i = 0; i < N1; i++ ) {
            // Get the baryon spinor in the B rest frame. Also just use u and not i*u,
            // since the imaginary constant factor is not needed for the probability
            EvtDiracSpinor u = baryon1->spParent( i );

            // Second baryon
            for ( int j = 0; j < N2; j++ ) {
                EvtDiracSpinor v = baryon2->spParent( j );

                // Hadronic currents
                std::vector<EvtVector4C> hadCurrents =
                    this->getHadronicCurrents( u, v, p, gMtmTerms, fMtmTerms );

                // First amplitude terms: 3rd current already has the form factor coeffs applied (gMtmTerms)
                EvtVector4C amp1 = g1 * sign1 * hadCurrents[0] +
                                   g2 * sign2 * hadCurrents[1] +
                                   sign3 * hadCurrents[2];

                // Second amplitude terms: 6th current already has the form factor coeffs applied (fMtmTerms)
                EvtVector4C amp2 = f1 * sign4 * hadCurrents[3] +
                                   f2 * sign5 * hadCurrents[4] +
                                   sign6 * hadCurrents[5];

                EvtVector4C hadAmp;
                if ( sameParity == 1 ) {
                    hadAmp = amp1 - amp2;
                } else {
                    hadAmp = amp2 - amp1;
                }

                amp.vertex( i, j, 0, l1 * hadAmp );
                amp.vertex( i, j, 1, l2 * hadAmp );

            }    // j

        }    // i

    } else if ( ( type1 == EvtSpinType::DIRAC &&
                  type2 == EvtSpinType::RARITASCHWINGER ) ||
                ( type1 == EvtSpinType::RARITASCHWINGER &&
                  type2 == EvtSpinType::DIRAC ) ) {
        // Handle the case of one Dirac-type daughter (not including the leptons), e.g. one proton, and one
        // Rarita-Schwinger-type (spin 3/2) daughter e.g. B -> p N(1520) l nu

        // Form factor parameters
        EvtBToDiBaryonlnupQCDFF::FormFactors FF;
        ffModel_.getRaritaFF( parent, m_dibaryon, FF );

        if ( sameParity == 1 ) {
            f1 = FF.F1;
            f2 = FF.F2;
            f3 = FF.F3;
            f4 = FF.F4;
            f5 = FF.F5;
            g1 = FF.G1;
            g2 = FF.G2;
            g3 = FF.G3;
            g4 = FF.G4;
            g5 = FF.G5;
        } else {
            // Swap coeffs: f_i <--> g_i
            f1 = FF.G1;
            f2 = FF.G2;
            f3 = FF.G3;
            f4 = FF.G4;
            f5 = FF.G5;
            g1 = FF.F1;
            g2 = FF.F2;
            g3 = FF.F3;
            g4 = FF.F4;
            g5 = FF.F5;
        }

        EvtVector4R gMtmTerms = g3 * p + g4 * pSum + g5 * pDiff;
        EvtVector4R fMtmTerms = f3 * p + f4 * pSum + f5 * pDiff;

        if ( type1 == EvtSpinType::DIRAC ) {
            // First baryon is Dirac
            for ( int i = 0; i < N1; i++ ) {
                // Get the baryon spinor in the B rest frame. Also just use u and not i*u,
                // since the imaginary constant factor is not needed for the probability
                EvtDiracSpinor u = baryon1->spParent( i );

                // Second baryon is RS-type
                for ( int j = 0; j < N2; j++ ) {
                    EvtRaritaSchwinger vRS = baryon2->spRSParent( j );

                    EvtDiracSpinor v;
                    for ( int k = 0; k < 4; k++ ) {
                        v.set_spinor( k, vRS.getVector( k ) * p0 );
                    }

                    // Hadronic currents
                    std::vector<EvtVector4C> hadCurrents =
                        this->getHadronicCurrents( u, v, p, gMtmTerms, fMtmTerms );

                    // First amplitude terms: 3rd current already has the form factor coeffs applied (gMtmTerms)
                    EvtVector4C amp1 = g1 * sign1 * hadCurrents[0] +
                                       g2 * sign2 * hadCurrents[1] +
                                       sign3 * hadCurrents[2];

                    // Second amplitude terms: 6th current already has the form factor coeffs applied (fMtmTerms)
                    EvtVector4C amp2 = f1 * sign4 * hadCurrents[3] +
                                       f2 * sign5 * hadCurrents[4] +
                                       sign6 * hadCurrents[5];

                    EvtVector4C hadAmp;
                    if ( sameParity == 1 ) {
                        hadAmp = amp1 - amp2;
                    } else {
                        hadAmp = amp2 - amp1;
                    }

                    amp.vertex( i, j, 0, l1 * hadAmp );
                    amp.vertex( i, j, 1, l2 * hadAmp );

                }    // j

            }    // i

        } else if ( type2 == EvtSpinType::DIRAC ) {
            // Same as before, but where the first daughter is RS-type, e.g. B -> N(1520) p l nu

            // First baryon is RS
            for ( int i = 0; i < N1; i++ ) {
                // Get the baryon spinor in the B rest frame
                EvtRaritaSchwinger uRS = baryon1->spRSParent( i );

                EvtDiracSpinor u;
                for ( int k = 0; k < 4; k++ ) {
                    u.set_spinor( k, uRS.getVector( k ) * p0 );
                }

                // Second baryon is Dirac
                for ( int j = 0; j < N2; j++ ) {
                    EvtDiracSpinor v = baryon2->spParent( j );

                    // Hadronic currents
                    std::vector<EvtVector4C> hadCurrents =
                        this->getHadronicCurrents( u, v, p, gMtmTerms, fMtmTerms );

                    // First amplitude terms: 3rd current already has the form factor coeffs applied (gMtmTerms)
                    EvtVector4C amp1 = g1 * sign1 * hadCurrents[0] +
                                       g2 * sign2 * hadCurrents[1] +
                                       sign3 * hadCurrents[2];

                    // Second amplitude terms: 6th current already has the form factor coeffs applied (fMtmTerms)
                    EvtVector4C amp2 = f1 * sign4 * hadCurrents[3] +
                                       f2 * sign5 * hadCurrents[4] +
                                       sign6 * hadCurrents[5];

                    EvtVector4C hadAmp;
                    if ( sameParity == 1 ) {
                        hadAmp = amp1 - amp2;
                    } else {
                        hadAmp = amp2 - amp1;
                    }

                    amp.vertex( i, j, 0, l1 * hadAmp );
                    amp.vertex( i, j, 1, l2 * hadAmp );

                }    // j

            }    // i

        }    // RS daughter check

    }    // Have Dirac and RS baryons
}

std::vector<EvtVector4C> EvtSLDiBaryonAmp::getHadronicCurrents(
    const EvtDiracSpinor& u, const EvtDiracSpinor& v, const EvtVector4R& p,
    const EvtVector4R& gMtmTerms, const EvtVector4R& fMtmTerms ) const
{
    // Store the currents used in Eq 6 (in order of appearance)
    std::vector<EvtVector4C> currents;
    currents.reserve( 6 );

    EvtDiracSpinor g5v = EvtGammaMatrix::g5() * v;

    // ubar*gamma*gamma5*v
    EvtVector4C current1 = EvtLeptonACurrent( u, v );
    currents.push_back( current1 );

    // ubar*sigma*p*gamma5*v -> [ubar*sigma*(gamma5*v)]*p
    EvtTensor4C TC1 = EvtLeptonTCurrent( u, g5v );
    // Contract tensor with 4-momentum
    EvtVector4C current2 = TC1.cont2( p );
    currents.push_back( current2 );

    // ubar*p*gamma5*v; "p" = p, pSum and pDiff
    EvtComplex PC1 = EvtLeptonPCurrent( u, v );
    EvtVector4C current3 = PC1 * gMtmTerms;
    currents.push_back( current3 );

    // ubar*gamma*v
    EvtVector4C current4 = EvtLeptonVCurrent( u, v );
    currents.push_back( current4 );

    // ubar*sigma*p*v -> [ubar*sigma*v]*p
    EvtTensor4C TC2 = EvtLeptonTCurrent( u, v );
    // Contract tensor with 4-momentum
    EvtVector4C current5 = TC2.cont2( p );
    currents.push_back( current5 );

    // ubar*p*v; "p" = p, pSum and pDiff
    EvtComplex S1 = EvtLeptonSCurrent( u, v );
    EvtVector4C current6 = S1 * fMtmTerms;
    currents.push_back( current6 );

    return currents;
}

int EvtSLDiBaryonAmp::checkDibaryonParity( const EvtId& id1, const EvtId& id2,
                                           const int J1, const int J2 ) const
{
    // Get intrisic parities of the two baryons, then multiply by (-1)^|J1 - J2|.
    // Note here that the J1 and J2 function arguments = 2*spin
    int par1 = this->getBaryonParity( id1 );
    int par2 = this->getBaryonParity( id2 );

    // mult should be either 0 or 1 for allowed Dirac/RS baryon pairs
    int mult = static_cast<int>( pow( -1.0, 0.5 * fabs( J1 - J2 ) ) );

    int dbParity = par1 * par2 * mult;

    // Initialise result to 1, i.e. dibaryon parity = B parity = -1
    int result( 1 );

    // Dibaryon parity is opposite to the negative B parity
    if ( dbParity > 0 ) {
        result = -1;
    }

    return result;
}

int EvtSLDiBaryonAmp::getBaryonParity( const EvtId& id ) const
{
    // Initialise parity to +1
    int parity( 1 );

    // List of baryons with parity = +1
    static EvtIdSet posParity( "p+", "Delta+", "Lambda_c+",
                               "anti-Lambda_c(2593)-", "anti-Lambda_c(2625)-",
                               "N(1440)+", "anti-N(1520)-", "anti-N(1535)-",
                               "anti-N(1650)-", "anti-N(1700)-", "N(1710)+",
                               "N(1720)+" );

    // If the baryon id is not in the list, set the parity to -1
    if ( !posParity.contains( id ) ) {
        parity = -1;
    }

    return parity;
}
