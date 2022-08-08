
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

#include "EvtGenBase/EvtGenKine.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <cmath>
#include <iostream>

using std::endl;

double EvtPawt( double a, double b, double c )
{
    double temp = ( a * a - ( b + c ) * ( b + c ) ) *
                  ( a * a - ( b - c ) * ( b - c ) );

    if ( temp <= 0 ) {
        return 0.0;
    }

    return sqrt( temp ) / ( 2.0 * a );
}

double EvtGenKine::PhaseSpace( int ndaug, double mass[30], EvtVector4R p4[30],
                               double mp )

//  N body phase space routine.  Send parent with
//  daughters already defined ( Number and masses )
//  Returns four vectors in parent frame.

{
    double energy, p3, alpha, beta;

    if ( ndaug == 1 ) {
        p4[0].set( mass[0], 0.0, 0.0, 0.0 );
        return 1.0;
    }

    if ( ndaug == 2 ) {
        //Two body phase space

        energy = ( mp * mp + mass[0] * mass[0] - mass[1] * mass[1] ) /
                 ( 2.0 * mp );

        p3 = 0.0;
        if ( energy > mass[0] ) {
            p3 = sqrt( energy * energy - mass[0] * mass[0] );
        }

        p4[0].set( energy, 0.0, 0.0, p3 );

        energy = mp - energy;
        p3 = -1.0 * p3;
        p4[1].set( energy, 0.0, 0.0, p3 );

        //Now rotate four vectors.

        alpha = EvtRandom::Flat( EvtConst::twoPi );
        beta = acos( EvtRandom::Flat( -1.0, 1.0 ) );

        p4[0].applyRotateEuler( alpha, beta, -alpha );
        p4[1].applyRotateEuler( alpha, beta, -alpha );

        return 1.0;
    }

    if ( ndaug != 2 ) {
        double wtmax = 0.0;
        double pm[5][30], pmin, pmax, psum, rnd[30];
        double ran, wt, pa, costh, sinth, phi, p[4][30], be[4], bep, temp;
        int i, il, ilr, i1, il1u, il1, il2r, ilu;
        int il2 = 0;

        for ( i = 0; i < ndaug; i++ ) {
            pm[4][i] = 0.0;
            rnd[i] = 0.0;
        }

        pm[0][0] = mp;
        pm[1][0] = 0.0;
        pm[2][0] = 0.0;
        pm[3][0] = 0.0;
        pm[4][0] = mp;

        psum = 0.0;
        for ( i = 1; i < ndaug + 1; i++ ) {
            psum = psum + mass[i - 1];
        }

        pm[4][ndaug - 1] = mass[ndaug - 1];

        switch ( ndaug ) {
            case 1:
                wtmax = 1.0 / 16.0;
                break;
            case 2:
                wtmax = 1.0 / 150.0;
                break;
            case 3:
                wtmax = 1.0 / 2.0;
                break;
            case 4:
                wtmax = 1.0 / 5.0;
                break;
            case 5:
                wtmax = 1.0 / 15.0;
                break;
            case 6:
                wtmax = 1.0 / 15.0;
                break;
            case 7:
                wtmax = 1.0 / 15.0;
                break;
            case 8:
                wtmax = 1.0 / 15.0;
                break;
            case 9:
                wtmax = 1.0 / 15.0;
                break;
            case 10:
                wtmax = 1.0 / 15.0;
                break;
            case 11:
                wtmax = 1.0 / 15.0;
                break;
            case 12:
                wtmax = 1.0 / 15.0;
                break;
            case 13:
                wtmax = 1.0 / 15.0;
                break;
            case 14:
                wtmax = 1.0 / 15.0;
                break;
            case 15:
                wtmax = 1.0 / 15.0;
                break;
            default:
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "too many daughters for phase space..." << ndaug << " "
                    << mp << endl;
                ;
                break;
        }

        pmax = mp - psum + mass[ndaug - 1];

        pmin = 0.0;

        for ( ilr = 2; ilr < ndaug + 1; ilr++ ) {
            il = ndaug + 1 - ilr;
            pmax = pmax + mass[il - 1];
            pmin = pmin + mass[il + 1 - 1];
            wtmax = wtmax * EvtPawt( pmax, pmin, mass[il - 1] );
        }

        do {
            rnd[0] = 1.0;
            il1u = ndaug - 1;

            for ( il1 = 2; il1 < il1u + 1; il1++ ) {
                ran = EvtRandom::Flat();
                for ( il2r = 2; il2r < il1 + 1; il2r++ ) {
                    il2 = il1 + 1 - il2r;
                    if ( ran <= rnd[il2 - 1] )
                        goto two39;
                    rnd[il2 + 1 - 1] = rnd[il2 - 1];
                }
            two39:
                rnd[il2 + 1 - 1] = ran;
            }

            rnd[ndaug - 1] = 0.0;
            wt = 1.0;
            for ( ilr = 2; ilr < ndaug + 1; ilr++ ) {
                il = ndaug + 1 - ilr;
                pm[4][il - 1] = pm[4][il + 1 - 1] + mass[il - 1] +
                                ( rnd[il - 1] - rnd[il + 1 - 1] ) * ( mp - psum );
                wt = wt *
                     EvtPawt( pm[4][il - 1], pm[4][il + 1 - 1], mass[il - 1] );
            }
            if ( wt > wtmax ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "wtmax to small in EvtPhaseSpace with " << ndaug
                    << " daughters" << endl;
                ;
            }
        } while ( wt < EvtRandom::Flat( wtmax ) );

        ilu = ndaug - 1;

        for ( il = 1; il < ilu + 1; il++ ) {
            pa = EvtPawt( pm[4][il - 1], pm[4][il + 1 - 1], mass[il - 1] );
            costh = EvtRandom::Flat( -1.0, 1.0 );
            sinth = sqrt( 1.0 - costh * costh );
            phi = EvtRandom::Flat( EvtConst::twoPi );
            p[1][il - 1] = pa * sinth * cos( phi );
            p[2][il - 1] = pa * sinth * sin( phi );
            p[3][il - 1] = pa * costh;
            pm[1][il + 1 - 1] = -p[1][il - 1];
            pm[2][il + 1 - 1] = -p[2][il - 1];
            pm[3][il + 1 - 1] = -p[3][il - 1];
            p[0][il - 1] = sqrt( pa * pa + mass[il - 1] * mass[il - 1] );
            pm[0][il + 1 - 1] = sqrt( pa * pa +
                                      pm[4][il + 1 - 1] * pm[4][il + 1 - 1] );
        }

        p[0][ndaug - 1] = pm[0][ndaug - 1];
        p[1][ndaug - 1] = pm[1][ndaug - 1];
        p[2][ndaug - 1] = pm[2][ndaug - 1];
        p[3][ndaug - 1] = pm[3][ndaug - 1];

        for ( ilr = 2; ilr < ndaug + 1; ilr++ ) {
            il = ndaug + 1 - ilr;
            be[0] = pm[0][il - 1] / pm[4][il - 1];
            be[1] = pm[1][il - 1] / pm[4][il - 1];
            be[2] = pm[2][il - 1] / pm[4][il - 1];
            be[3] = pm[3][il - 1] / pm[4][il - 1];

            for ( i1 = il; i1 < ndaug + 1; i1++ ) {
                bep = be[1] * p[1][i1 - 1] + be[2] * p[2][i1 - 1] +
                      be[3] * p[3][i1 - 1] + be[0] * p[0][i1 - 1];
                temp = ( p[0][i1 - 1] + bep ) / ( be[0] + 1.0 );
                p[1][i1 - 1] = p[1][i1 - 1] + temp * be[1];
                p[2][i1 - 1] = p[2][i1 - 1] + temp * be[2];
                p[3][i1 - 1] = p[3][i1 - 1] + temp * be[3];
                p[0][i1 - 1] = bep;
            }
        }

        for ( ilr = 0; ilr < ndaug; ilr++ ) {
            p4[ilr].set( p[0][ilr], p[1][ilr], p[2][ilr], p[3][ilr] );
        }

        return 1.0;
    }

    return 1.0;
}

double EvtGenKine::PhaseSpacePole( double M, double m1, double m2, double m3,
                                   double a, EvtVector4R p4[10] )

//  generate kinematics for 3 body decays, pole for the m1,m2 mass.

{
    //f1   = 1  (phasespace)
    //f2   = a*(1/m12sq)^2

    double m12sqmax = ( M - m3 ) * ( M - m3 );
    double m12sqmin = ( m1 + m2 ) * ( m1 + m2 );

    double m13sqmax = ( M - m2 ) * ( M - m2 );
    double m13sqmin = ( m1 + m3 ) * ( m1 + m3 );

    double v1 = ( m12sqmax - m12sqmin ) * ( m13sqmax - m13sqmin );
    double v2 = a * ( 1.0 / m12sqmin - 1.0 / m12sqmax ) * ( m13sqmax - m13sqmin );

    double m12sq, m13sq;

    double r = v1 / ( v1 + v2 );

    double m13min, m13max;

    do {
        m13sq = EvtRandom::Flat( m13sqmin, m13sqmax );

        if ( r > EvtRandom::Flat() ) {
            m12sq = EvtRandom::Flat( m12sqmin, m12sqmax );
        } else {
            m12sq = 1.0 /
                    ( 1.0 / m12sqmin -
                      EvtRandom::Flat() * ( 1.0 / m12sqmin - 1.0 / m12sqmax ) );
        }

        //kinematically allowed?
        double E3star = ( M * M - m12sq - m3 * m3 ) / sqrt( 4 * m12sq );
        double E1star = ( m12sq + m1 * m1 - m2 * m2 ) / sqrt( 4 * m12sq );
        double p3star = sqrt( E3star * E3star - m3 * m3 );
        double p1star = sqrt( E1star * E1star - m1 * m1 );
        m13max = ( E3star + E1star ) * ( E3star + E1star ) -
                 ( p3star - p1star ) * ( p3star - p1star );
        m13min = ( E3star + E1star ) * ( E3star + E1star ) -
                 ( p3star + p1star ) * ( p3star + p1star );

    } while ( m13sq < m13min || m13sq > m13max );

    double E2 = ( M * M + m2 * m2 - m13sq ) / ( 2.0 * M );
    double E3 = ( M * M + m3 * m3 - m12sq ) / ( 2.0 * M );
    double E1 = M - E2 - E3;
    double p1mom = sqrt( E1 * E1 - m1 * m1 );
    double p3mom = sqrt( E3 * E3 - m3 * m3 );
    double cost13 = ( 2.0 * E1 * E3 + m1 * m1 + m3 * m3 - m13sq ) /
                    ( 2.0 * p1mom * p3mom );

    //EvtGenReport(EVTGEN_INFO,"EvtGen") << m13sq << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << m12sq << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << E1 << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << E2 << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << E3 << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << p1mom << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << p3mom << endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << cost13 << endl;

    p4[2].set( E3, 0.0, 0.0, p3mom );
    p4[0].set( E1, p1mom * sqrt( 1.0 - cost13 * cost13 ), 0.0, p1mom * cost13 );
    p4[1].set( E2, -p1mom * sqrt( 1.0 - cost13 * cost13 ), 0.0,
               -p1mom * cost13 - p3mom );

    //EvtGenReport(EVTGEN_INFO,"EvtGen") << "p4:"<<p4[0]<<p4[1]<<p4[2]<<endl;

    double alpha = EvtRandom::Flat( EvtConst::twoPi );
    double beta = acos( EvtRandom::Flat( -1.0, 1.0 ) );
    double gamma = EvtRandom::Flat( EvtConst::twoPi );

    p4[0].applyRotateEuler( alpha, beta, gamma );
    p4[1].applyRotateEuler( alpha, beta, gamma );
    p4[2].applyRotateEuler( alpha, beta, gamma );

    return 1.0 + a / ( m12sq * m12sq );
}

/*
 * Function which takes two invariant masses squared in 3-body decay and
 * parent after makeDaughters() and generateMassTree() and
 * calculates/generates momenta of daughters and sets those.
 */
void EvtGenKine::ThreeBodyKine( const double m12Sq, const double m23Sq,
                                EvtParticle* p )
{
    const double mParent = p->mass();
    EvtParticle* daug1 = p->getDaug( 0 );
    EvtParticle* daug2 = p->getDaug( 1 );
    EvtParticle* daug3 = p->getDaug( 2 );
    const double mDaug1 = daug1->mass();
    const double mDaug2 = daug2->mass();
    const double mDaug3 = daug3->mass();
    const double mParentSq{ mParent * mParent };
    const double mDaug1Sq{ mDaug1 * mDaug1 };
    const double mDaug2Sq{ mDaug2 * mDaug2 };
    const double mDaug3Sq{ mDaug3 * mDaug3 };
    const double invMParent{ 1. / mParent };

    const double En1 = 0.5 * ( mParentSq + mDaug1Sq - m23Sq ) * invMParent;
    const double En3 = 0.5 * ( mParentSq + mDaug3Sq - m12Sq ) * invMParent;
    const double En2 = mParent - En1 - En3;
    const double p1mag = std::sqrt( En1 * En1 - mDaug1Sq );
    const double p2mag = std::sqrt( En2 * En2 - mDaug2Sq );
    double cosPhi = 0.5 * ( mDaug1Sq + mDaug2Sq + 2 * En1 * En2 - m12Sq ) /
                    ( p1mag * p2mag );

    double sinPhi = std::sqrt( 1 - cosPhi * cosPhi );
    if ( EvtRandom::Flat( 0., 1. ) > 0.5 ) {
        sinPhi *= -1;
    }
    const double p2x = p2mag * cosPhi;
    const double p2y = p2mag * sinPhi;
    const double p3x = -p1mag - p2x;
    const double p3y = -p2y;

    // Construct 4-momenta and rotate them randomly in space
    EvtVector4R p1( En1, p1mag, 0., 0. );
    EvtVector4R p2( En2, p2x, p2y, 0. );
    EvtVector4R p3( En3, p3x, p3y, 0. );
    const double euler1 = EvtRandom::Flat( 0., EvtConst::twoPi );
    const double euler2 = std::acos( EvtRandom::Flat( -1.0, 1.0 ) );
    const double euler3 = EvtRandom::Flat( 0., EvtConst::twoPi );
    p1.applyRotateEuler( euler1, euler2, euler3 );
    p2.applyRotateEuler( euler1, euler2, euler3 );
    p3.applyRotateEuler( euler1, euler2, euler3 );

    // set momenta for daughters
    daug1->init( daug1->getId(), p1 );
    daug2->init( daug2->getId(), p2 );
    daug3->init( daug3->getId(), p3 );

    return;
}
