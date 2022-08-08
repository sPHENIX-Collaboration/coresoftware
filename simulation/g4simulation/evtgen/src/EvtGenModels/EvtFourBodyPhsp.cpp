
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

#include <cmath>

#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtFourBodyPhsp.hh"

std::string EvtFourBodyPhsp::getName()
{
    return "FOURBODYPHSP";
}

EvtDecayBase* EvtFourBodyPhsp::clone()
{
    return new EvtFourBodyPhsp;
}

std::array<double, 4> EvtFourBodyPhsp::phspFactor(
    const double mM, const double m12, const double m34,
    std::array<double, 4>& daughters ) const
{
    std::array<double, 4> result;
    result[1] = twoBodyMomentum( mM, m12, m34 );
    result[2] = twoBodyMomentum( m12, daughters[0], daughters[1] );
    result[3] = twoBodyMomentum( m34, daughters[2], daughters[3] );
    result[0] = result[1] * result[2] * result[3];

    return result;
}

void EvtFourBodyPhsp::init()
{
    // Check that we have right number of daughters
    checkNDaug( 4 );

    // Check whether mother is quasi-stable
    auto parent = getParentId();
    double width = EvtPDL::getWidth( parent );
    if ( width > 1e-6 ) {
        m_stableMother = false;
    }

    // Check whether all daughters are stable
    for ( int i = 0; i < 4; ++i ) {
        auto daughter = getDaug( i );
        width = EvtPDL::getWidth( daughter );
        if ( width > 1e-6 ) {
            m_stableDaughters = false;
            m_daughterMasses[i] = EvtPDL::getMinMass( daughter );
        } else {
            m_daughterMasses[i] = EvtPDL::getMass( daughter );
        }
    }

    // check correct number of arguments
    checkNArg( 0, 2, 4 );
    double mass1 = m_daughterMasses[0];
    double mass2 = m_daughterMasses[1];
    double mass3 = m_daughterMasses[2];
    double mass4 = m_daughterMasses[3];
    double mMother = EvtPDL::getMaxMass( parent );
    if ( getNArg() > 2 ) {
        m_m12Min = getArg( 0 );
        m_m12Max = getArg( 1 );
        m_m34Min = getArg( 2 );
        m_m34Max = getArg( 3 );
    } else {
        if ( getNArg() > 0 ) {
            m_m12Min = getArg( 0 );
            m_m12Max = getArg( 1 );
        } else {
            m_m12Min = mass1 + mass2;
            m_m12Max = mMother - mass3 - mass4;
        }
        m_m34Min = mass3 + mass4;
        m_m34Max = mMother - mass1 - mass2;
        if ( m_stableDaughters == false || m_stableMother == false ) {
            m_fixedBoundary = false;
        }
    }
    // Make sure that we have correct boundaries
    if ( m_m12Min < mass1 + mass2 ) {
        m_m12Min = mass1 + mass2;
    }
    if ( m_m12Max > mMother - mass3 - mass4 ) {
        m_m12Max = mMother - mass3 - mass4;
    }
    if ( m_m34Min < mass3 + mass4 ) {
        m_m34Min = mass3 + mass4;
    }
    if ( m_m34Max > mMother - mass1 - mass2 ) {
        m_m34Max = mMother - mass1 - mass2;
    }

    if ( m_stableDaughters && m_stableMother ) {
        m_boundaryShape = determineBoundaryShape( m_m12Min, m_m12Max, m_m34Max,
                                                  mMother );
    } else {
        m_boundaryShape = Shape::variable;
    }
    // If we have fixed boundary, we can precalculate some variables for
    // m12 and m34 generation
    if ( m_fixedBoundary ) {
        if ( m_boundaryShape == Shape::trapezoid ) {
            const double m12Diff = m_m12Max - m_m12Min;
            const double minSum = m_m12Min + m_m34Min;
            m_trapNorm = ( mMother - m_m34Min ) * m12Diff -
                         0.5 * ( m12Diff * ( m_m12Max + m_m12Min ) );
            m_trapCoeff1 = mMother - m_m34Min;
            m_trapCoeff2 = mMother * mMother - 2 * mMother * minSum +
                           minSum * minSum;
        }
        if ( m_boundaryShape == Shape::pentagon ) {
            m_pentagonSplit = mMother - m_m34Max;
            const double area1 = ( m_pentagonSplit - m_m12Min ) *
                                 ( m_m34Max - m_m34Min );
            const double pm12Diff = m_m12Max - m_pentagonSplit;
            const double area2 = 0.5 * pm12Diff *
                                     ( mMother + m_m34Max - m_m12Max ) -
                                 pm12Diff * m_m34Min;
            m_pentagonFraction = area1 / ( area1 + area2 );
            const double m12Diff = m_m12Max - m_pentagonSplit;
            const double minSum = m_pentagonSplit + m_m34Min;
            m_trapNorm = ( mMother - m_m34Min ) * m12Diff -
                         0.5 * ( m12Diff * ( m_m12Max + m_pentagonSplit ) );
            m_trapCoeff1 = mMother - m_m34Min;
            m_trapCoeff2 = mMother * mMother - 2 * mMother * minSum +
                           minSum * minSum;
        }
    }
}

void EvtFourBodyPhsp::initProbMax()
{
    double startM12 = m_m12Min + ( m_m12Max - m_m12Min ) / 20.;
    double startM34 = m_m34Min + ( m_m34Max - m_m34Min ) / 20.;
    bool contCond = true;
    int iteration = 0;

    auto parent = getParentId();
    double mMother = EvtPDL::getMaxMass( parent );

    double funcValue = 0;
    while (contCond){
        ++iteration;
        double currentM12 = startM12;
        double currentM34 = startM34;
        funcValue = phspFactor( mMother, currentM12, currentM34,
                                m_daughterMasses )[0];
        // Maximum along m12
        double step = ( m_m12Max - m_m12Min ) / 100.;
        while ( step > 1e-4 ) {
            double point1 = currentM12 + step;
            if ( point1 > m_m12Max ) {
                point1 = m_m12Max;
            }
            if ( currentM34 > mMother - point1 ) {
                point1 = mMother - currentM34;
            }
            double point2 = currentM12 - step;
            if ( point2 < m_m12Min ) {
                point2 = m_m12Min;
            }
            double value1 = phspFactor( mMother, point1, currentM34,
                                        m_daughterMasses )[0];
            double value2 = phspFactor( mMother, point2, currentM34,
                                        m_daughterMasses )[0];
            if ( value1 > funcValue && value1 > value2 ) {
                currentM12 = point1;
                funcValue = value1;
            } else if ( value2 > funcValue ) {
                currentM12 = point2;
                funcValue = value2;
            }
            step /= 2.;
        }
        // Maximum along m34
        step = ( mMother - currentM12 - m_m34Min ) / 100.;
        while ( step > 1e-4 ) {
            double point1 = currentM34 + step;
            if ( point1 > m_m34Max ) {
                point1 = m_m34Max;
            }
            if ( point1 > mMother - currentM12 ) {
                point1 = mMother - currentM12;
            }
            double point2 = currentM34 - step;
            if ( point2 < m_m34Min ) {
                point2 = m_m34Min;
            }
            double value1 = phspFactor( mMother, currentM12, point1,
                                        m_daughterMasses )[0];
            double value2 = phspFactor( mMother, currentM12, point2,
                                        m_daughterMasses )[0];
            if ( value1 > funcValue && value1 > value2 ) {
                currentM34 = point1;
                funcValue = value1;
            } else if ( value2 > funcValue ) {
                currentM34 = point2;
                funcValue = value2;
            }
            step /= 2.;
        }

        // Check termination condition
       double m12Diff = currentM12 - startM12;
       double m34Diff = currentM34 - startM34;
       double distSq = m12Diff * m12Diff + m34Diff * m34Diff;
       if (distSq < 1e-8 || iteration > 50){
           contCond = false;
       }
       startM12 = currentM12;
       startM34 = currentM34;
    }

    setProbMax( funcValue * 1.05 );
}

void EvtFourBodyPhsp::decay( EvtParticle* parent )
{

    parent->makeDaughters( getNDaug(), getDaugs() );
    bool massTreeStatus = parent->generateMassTree();
    if ( !massTreeStatus ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Failed to generate daughters masses in EvtFourBodyPhsp."
            << std::endl;
        ::abort();
    }

    double mMother = parent->mass();

    // Need to check whether boundaries are OK and whether we need to work
    // out boundary shape
    double cM12Min, cM12Max;
    double cM34Min, cM34Max;
    EvtFourBodyPhsp::Shape cShape;
    if ( m_fixedBoundary ) {
        cM12Min = m_m12Min;
        cM12Max = m_m12Max;
        cM34Min = m_m34Min;
        cM34Max = m_m34Max;
        cShape = m_boundaryShape;
    } else {
        // In this case at least one particle has non-zero width and thus
        // boundaries and shape of the region can change
        for ( int i = 0; i < 4; ++i ) {
            auto daughter = parent->getDaug( i );
            m_daughterMasses[i] = daughter->mass();
        }
        cM12Min = m_m12Min > ( m_daughterMasses[0] + m_daughterMasses[1] )
                      ? m_m12Min
                      : m_daughterMasses[0] + m_daughterMasses[1];
        cM12Max = m_m12Max <
                          ( mMother - m_daughterMasses[2] - m_daughterMasses[3] )
                      ? m_m12Max
                      : mMother - m_daughterMasses[2] - m_daughterMasses[3];
        cM34Min = m_m34Min > ( m_daughterMasses[2] + m_daughterMasses[3] )
                      ? m_m34Min
                      : m_daughterMasses[2] + m_daughterMasses[3];
        cM34Max = m_m34Max <
                          ( mMother - m_daughterMasses[0] - m_daughterMasses[1] )
                      ? m_m34Max
                      : mMother - m_daughterMasses[0] - m_daughterMasses[1];
        cShape = determineBoundaryShape( cM12Min, cM12Max, cM34Max, mMother );
    }

    // Generate m12 and m34
    auto masses = generatePairMasses( cM12Min, cM12Max, cM34Min, cM34Max,
                                      mMother, cShape );
    const double m12 = masses.first;
    const double m34 = masses.second;

    // calculate probability, it will return array with 4 elements with
    // probability, q, p1 and p3
    auto probEval = phspFactor( mMother, m12, m34, m_daughterMasses );
    setProb( probEval[0] );

    // initialise kinematics
    const double cosTheta1 = EvtRandom::Flat(-1.0, 1.0);
    const double sinTheta1 = std::sqrt( 1 - cosTheta1 * cosTheta1 );
    const double cosTheta3 = EvtRandom::Flat(-1.0, 1.0);
    const double sinTheta3 = std::sqrt( 1 - cosTheta3 * cosTheta3 );
    const double phi = EvtRandom::Flat( 0., EvtConst::twoPi );
    // m12 and m34 are put along z-axis, 1 and 2 go to x-z plane and 3-4
    // plane is rotated by phi compared to 1-2 plane. All momenta are set
    // in 12 and 34 rest frames and then boosted to parent rest frame
    const double p1x = probEval[2] * sinTheta1;
    const double p1z = probEval[2] * cosTheta1;
    const double p1Sq = probEval[2] * probEval[2];
    const double en1 = std::sqrt( m_daughterMasses[0] * m_daughterMasses[0] +
                                  p1Sq );
    const double en2 = std::sqrt( m_daughterMasses[1] * m_daughterMasses[1] +
                                  p1Sq );
    const double p3T = probEval[3] * sinTheta3;
    const double p3x = p3T * std::cos( phi );
    const double p3y = p3T * std::sin( phi );
    const double p3z = probEval[3] * cosTheta3;
    const double p3Sq = probEval[3] * probEval[3];
    const double en3 = std::sqrt( m_daughterMasses[2] * m_daughterMasses[2] +
                                  p3Sq );
    const double en4 = std::sqrt( m_daughterMasses[3] * m_daughterMasses[3] +
                                  p3Sq );

    EvtVector4R mom1( en1, p1x, 0.0, p1z );
    EvtVector4R mom2( en2, -p1x, 0.0, -p1z );
    EvtVector4R mom3( en3, p3x, p3y, p3z );
    EvtVector4R mom4( en4, -p3x, -p3y, -p3z );

    const double qSq = probEval[1] * probEval[1];
    const double en12 = std::sqrt( m12 * m12 + qSq );
    const double en34 = std::sqrt( m34 * m34 + qSq );
    EvtVector4R q12( en12, 0.0, 0.0, probEval[1] );
    EvtVector4R q34( en34, 0.0, 0.0, -probEval[1] );
    mom1.applyBoostTo( q12 );
    mom2.applyBoostTo( q12 );
    mom3.applyBoostTo( q34 );
    mom4.applyBoostTo( q34 );

    // As final step, rotate everything randomly in space
    const double euler1 = EvtRandom::Flat( 0., EvtConst::twoPi );
    const double euler2 = std::acos( EvtRandom::Flat( -1.0, 1.0 ) );
    const double euler3 = EvtRandom::Flat( 0., EvtConst::twoPi );
    mom1.applyRotateEuler( euler1, euler2, euler3 );
    mom2.applyRotateEuler( euler1, euler2, euler3 );
    mom3.applyRotateEuler( euler1, euler2, euler3 );
    mom4.applyRotateEuler( euler1, euler2, euler3 );

    // Set momenta for daughters
    auto daug = parent->getDaug( 0 );
    daug->init( daug->getId(), mom1 );
    daug = parent->getDaug( 1 );
    daug->init( daug->getId(), mom2 );
    daug = parent->getDaug( 2 );
    daug->init( daug->getId(), mom3 );
    daug = parent->getDaug( 3 );
    daug->init( daug->getId(), mom4 );
}

EvtFourBodyPhsp::Shape EvtFourBodyPhsp::determineBoundaryShape(
    const double m12Min, const double m12Max, const double m34Max,
    const double mMother ) const
{
    double maxY = mMother - m12Min;
    const bool corner1 = m34Max < maxY;
    maxY = mMother - m12Max;
    const bool corner2 = m34Max < maxY;

    if ( corner1 && corner2 ) {
        return Shape::rectangle;
    } else if ( !corner1 && !corner2 ) {
        return Shape::trapezoid;
    }
    return Shape::pentagon;
}

std::pair<double, double> EvtFourBodyPhsp::generatePairMasses(
    const double m12Min, const double m12Max, const double m34Min,
    const double m34Max, const double mMother,
    const EvtFourBodyPhsp::Shape shape ) const
{
    switch ( shape ) {
        case EvtFourBodyPhsp::Shape::rectangle:
            return generateRectangle( m12Min, m12Max, m34Min, m34Max );
            break;
        case EvtFourBodyPhsp::Shape::trapezoid:
            return generateTrapezoid( m12Min, m12Max, m34Min, mMother );
            break;
        case EvtFourBodyPhsp::Shape::pentagon:
            double split, fraction;
            if ( m_fixedBoundary ) {
                split = m_pentagonSplit;
                fraction = m_pentagonFraction;
            } else {
                split = mMother - m34Max;
                const double area1 = ( split - m12Min ) * ( m34Max - m34Min );
                const double pm12Diff = m12Max - split;
                const double area2 = 0.5 * pm12Diff *
                                         ( mMother + m34Max - m12Max ) -
                                     pm12Diff * m34Min;
                fraction = area1 / ( area1 + area2 );
            }
            if ( EvtRandom::Flat() < fraction ) {
                return generateRectangle( m12Min, split, m34Min, m34Max );
            } else {
                return generateTrapezoid( split, m12Max, m34Min, mMother );
            }
            break;
        default:
            return std::make_pair( m12Min, m34Min );
            break;
    }
}

std::pair<double, double> EvtFourBodyPhsp::generateRectangle(
    const double m12Min, const double m12Max, const double m34Min,
    const double m34Max ) const
{
    return std::make_pair( EvtRandom::Flat( m12Min, m12Max ),
                           EvtRandom::Flat( m34Min, m34Max ) );
}

std::pair<double, double> EvtFourBodyPhsp::generateTrapezoid(
    const double m12Min, const double m12Max, const double m34Min,
    const double mMother ) const
{
    double norm, coeff1, coeff2;
    if ( m_fixedBoundary ) {
        norm = m_trapNorm;
        coeff1 = m_trapCoeff1;
        coeff2 = m_trapCoeff2;
    } else {
        const double m12Diff = m12Max - m12Min;
        const double minSum = m12Min + m34Min;
        norm = ( mMother - m34Min ) * m12Diff -
               0.5 * ( m12Diff * ( m12Max + m12Min ) );
        coeff1 = mMother - m34Min;
        coeff2 = mMother * mMother - 2 * mMother * minSum + minSum * minSum;
    }
    const double rnd = EvtRandom::Flat();
    const double m12 = coeff1 - std::sqrt( -2.0 * rnd * norm + coeff2 );
    const double m34 = EvtRandom::Flat( m34Min, mMother - m12 );
    return std::make_pair( m12, m34 );
}
