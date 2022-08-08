
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

#include "EvtGenModels/EvtHypNonLepton.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

EvtDecayBase* EvtHypNonLepton::clone()
{
    return new EvtHypNonLepton;
}

std::string EvtHypNonLepton::getName()
{
    return "HypNonLepton";
}

void EvtHypNonLepton::init()
{
    if ( getNArg() < 2 || getNArg() > 3 ) {    // alpha phi gamma delta
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton generator expected 2 or 3 arguments but found: "
            << getNArg() << std::endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "  1. Decay asymmetry parameter - alpha" << std::endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "  2. Parameter phi - in degrees (not radians)" << std::endl;
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "  3. Note on every x-th decay" << std::endl;
        ::abort();
    }

    if ( getNDaug() != 2 ) {    // Check that there are 2 daughters only
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton generator expected 2 daughters but found: "
            << getNDaug() << std::endl;
        ::abort();
    }

    // Check particles spins
    if ( EvtSpinType::getSpin2( EvtPDL::getSpinType( getParentId() ) ) != 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton generator expected dirac parent particle, but found "
            << EvtSpinType::getSpin2( EvtPDL::getSpinType( getParentId() ) )
            << " spin degrees of freedom" << std::endl;
        ::abort();
    }
    if ( EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 0 ) ) ) != 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton generator expected the first child to be dirac particle, but found "
            << EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 0 ) ) )
            << " spin degrees of freedom" << std::endl;
        ::abort();
    }
    if ( EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 1 ) ) ) != 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton generator expected the second child to be scalar particle, but found "
            << EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 1 ) ) )
            << " spin degrees of freedom" << std::endl;
        ::abort();
    }

    // Read all parameters
    m_alpha = getArg( 0 );
    m_phi = getArg( 1 ) * EvtConst::pi / 180;
    if ( getNArg() == 3 )
        m_noTries = static_cast<long>( getArg( 2 ) );
    else
        m_noTries = 0;

    // calculate additional parameters
    double p, M, m1, m2;
    double p_to_s, beta, delta, gamma;

    M = EvtPDL::getMass( getParentId() );
    m1 = EvtPDL::getMass( getDaug( 0 ) );
    m2 = EvtPDL::getMass( getDaug( 1 ) );

    if ( m1 + m2 >= M ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton found impossible decay: " << M
            << " --> " << m1 << " + " << m2 << " GeV\n"
            << std::endl;
        ::abort();
    }

    p = sqrt( M * M - ( m1 + m2 ) * ( m1 + m2 ) ) *
        sqrt( M * M - ( m1 - m2 ) * ( m1 - m2 ) ) / 2. / M;

    beta = sqrt( 1. - m_alpha * m_alpha ) * sin( m_phi );
    delta = -atan2( beta, m_alpha );
    gamma = sqrt( 1. - m_alpha * m_alpha - beta * beta );
    p_to_s = sqrt( ( 1. - gamma ) / ( 1. + gamma ) );

    m_B_to_A = p_to_s * ( m1 + sqrt( p * p + m1 * m1 ) ) / p *
               EvtComplex( cos( delta ), sin( delta ) );
}

void EvtHypNonLepton::initProbMax()
{
    double maxProb, m1, m2, M, p;

    M = EvtPDL::getMass( getParentId() );
    m1 = EvtPDL::getMass( getDaug( 0 ) );
    m2 = EvtPDL::getMass( getDaug( 1 ) );

    if ( m1 + m2 >= M ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << " ERROR: EvtHypNonLepton found impossible decay: " << M
            << " --> " << m1 << " + " << m2 << " GeV\n"
            << std::endl;
        ::abort();
    }

    p = sqrt( M * M - ( m1 + m2 ) * ( m1 + m2 ) ) *
        sqrt( M * M - ( m1 - m2 ) * ( m1 - m2 ) ) / 2 / M;
    maxProb = 16 * M *
              ( sqrt( p * p + m1 * m1 ) + m1 +
                abs( m_B_to_A ) * abs( m_B_to_A ) *
                    ( sqrt( p * p + m1 * m1 ) - m1 ) );
    //maxProb *= G_F*M_pi*M_pi;

    setProbMax( maxProb );
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << " EvtHypNonLepton set up maximum probability to " << maxProb
        << std::endl;
}

void EvtHypNonLepton::decay( EvtParticle* parent )
{
    parent->initializePhaseSpace( getNDaug(), getDaugs() );
    calcAmp( &_amp2, parent );
}

void EvtHypNonLepton::calcAmp( EvtAmp* amp, EvtParticle* parent )
{
    static long noTries = 0;
    int i;
    EvtComplex Matrix[2][2];

    //G_F  = 1.16637e-5;
    //M_pi = 0.13957;

    for ( i = 0; i < 4; i++ ) {
        //std::cout << "--------------------------------------------------" << std::endl;
        Matrix[i / 2][i % 2] = EvtLeptonSCurrent(
            parent->sp( i / 2 ), parent->getDaug( 0 )->spParent( i % 2 ) );
        //std::cout << "Matrix = " << Matrix[i/2][i%2] << std::endl;
        Matrix[i / 2][i % 2] -=
            m_B_to_A *
            EvtLeptonPCurrent( parent->sp( i / 2 ),
                               parent->getDaug( 0 )->spParent( i % 2 ) );
        //std::cout << "Matrix = " << Matrix[i/2][i%2] << std::endl;
        //Matrix[i/2][i%2] *= G_F*M_pi*M_pi;
        //std::cout << "Matrix = " << Matrix[i/2][i%2] << std::endl;
        //std::cout << "--------------------------------------------------" << std::endl;
        amp->vertex( i / 2, i % 2, Matrix[i / 2][i % 2] );
    }

    if ( m_noTries > 0 )
        if ( !( ( ++noTries ) % m_noTries ) )
            EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
                << " EvtHypNonLepton already finished " << noTries
                << " matrix element calculations" << std::endl;
}
