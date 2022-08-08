
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

#include "EvtGenModels/EvtBsquark.hh"

#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <string>

std::string EvtBsquark::getName()
{
    return "BSQUARK";
}

EvtDecayBase* EvtBsquark::clone()
{
    return new EvtBsquark;
}

void EvtBsquark::init()
{
    // check that there are 5 arguments
    checkNArg( 5 );
}

void EvtBsquark::initProbMax()
{
    //For now do not set a maximum.

    //SetProbMax(0.000000000005);
}

void EvtBsquark::decay( EvtParticle* p )
{
    static EvtId cquark = EvtPDL::getId( "c" );
    static EvtId anticquark = EvtPDL::getId( "anti-c" );

    static EvtIdSet leptons( "e-", "mu-", "tau-" );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    int charge = 1;

    EvtParticle* lepton;
    lepton = p->getDaug( 1 );
    if ( leptons.contains( lepton->getId() ) ) {
        charge = -1;
    }

    EvtDiracParticle charmquark;

    //this is a very crude approximation...
    if ( charge == -1 ) {
        charmquark.init( cquark, p->getDaug( 0 )->getP4() );
    } else {
        charmquark.init( anticquark, p->getDaug( 0 )->getP4() );
    }

    EvtVector4R p4c = p->getDaug( 0 )->getP4();

    EvtVector4R p4sn = p->getDaug( 2 )->getP4();

    EvtVector4R p4b( p->mass(), 0.0, 0.0, 0.0 );

    EvtComplex M[2][2];

    int il, ic;

    //project out the right handed current
    EvtGammaMatrix PR = 0.5 * ( EvtGammaMatrix::id() + EvtGammaMatrix::g5() );

    double tanbeta = getArg( 1 );
    double cosbeta = cos( atan( tanbeta ) );
    double sinbeta = sin( atan( tanbeta ) );

    double mb = 4.9;
    double mc = 1.3;
    double mw = 80.4;

    double Mass = getArg( 2 );
    double mu = getArg( 3 );
    double mchargino = getArg( 4 );

    double tan2phim = 2 * sqrt( 2.0 ) * mw * ( mu * cosbeta + Mass * sinbeta ) /
                      ( Mass * Mass - mu * mu +
                        2 * mw * mw * cos( 2 * atan( tanbeta ) ) );

    double phim = 0.5 * atan( tan2phim );

    EvtComplex U11 = cos( phim );
    EvtComplex U12 = sin( phim );
    EvtComplex U21 = -sin( phim );
    EvtComplex U22 = cos( phim );

    double tan2phip = 2 * sqrt( 2.0 ) * mw * ( mu * cosbeta + Mass * sinbeta ) /
                      ( Mass * Mass - mu * mu -
                        2 * mw * mw * cos( 2 * atan( tanbeta ) ) );

    double phip = 0.5 * atan( tan2phip );

    EvtComplex V11 = cos( phip );
    EvtComplex V12 = sin( phip );
    EvtComplex V21 = -sin( phip );
    EvtComplex V22 = cos( phip );

    double theta = getArg( 0 );
    double ctheta = cos( theta );
    double stheta = sin( theta );

    double vcsb = 0.08;
    double mchi1 = mchargino;
    double mchi2 = mchargino;

    //overall scale factor
    double g = 1.0;

    EvtComplex a1 = mchi1 * ( U11 * ctheta - mb * U12 * stheta /
                                                 ( sqrt( 2.0 ) * mw * cosbeta ) );
    EvtComplex a2 = mchi2 * ( U21 * ctheta - mb * U22 * stheta /
                                                 ( sqrt( 2.0 ) * mw * cosbeta ) );

    EvtComplex b1 = mc * conj( V12 ) * ctheta / ( sqrt( 2.0 ) * mw * sinbeta );
    EvtComplex b2 = mc * conj( V22 ) * ctheta / ( sqrt( 2.0 ) * mw * sinbeta );

    EvtComplex f1 = -( g * g * V11 * vcsb ) /
                    ( ( p4b - p4c ).mass2() - mchi1 * mchi1 );
    EvtComplex f2 = -( g * g * V21 * vcsb ) /
                    ( ( p4b - p4c ).mass2() - mchi1 * mchi2 );

    //EvtGenReport(EVTGEN_INFO,"EvtGen") <<g<<" "<<V11<<" "<<FL<<" "<<vcsb<<" "<<mchi1<<endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << "f1:"<<f1<<" "<<(p4b-p4c).mass2()<<endl;
    //EvtGenReport(EVTGEN_INFO,"EvtGen") << "f2:"<<f2<<" "<<(p4b-p4c).mass2()<<endl;

    //EvtGenReport(EVTGEN_INFO,"EvtGen") << "p4sn:"<<p4sn<<endl;

    EvtGammaMatrix pslash = p4sn.get( 0 ) * EvtGammaMatrix::g0() -
                            p4sn.get( 1 ) * EvtGammaMatrix::g1() -
                            p4sn.get( 2 ) * EvtGammaMatrix::g2() -
                            p4sn.get( 3 ) * EvtGammaMatrix::g3();

    //EvtGenReport(EVTGEN_INFO,"EvtGen") << "pslash:"<<pslash<<endl;

    for ( il = 0; il < 2; il++ ) {
        for ( ic = 0; ic < 2; ic++ ) {
            EvtComplex a = 0.0;
            EvtComplex b = 0.0;

            if ( charge == -1 ) {
                a = charmquark.spParent( ic ) * ( PR * lepton->spParent( il ) );
                b = charmquark.spParent( ic ) *
                    ( ( pslash * PR ) * lepton->spParent( il ) );
            } else {
                a = lepton->spParent( il ) * ( PR * charmquark.spParent( ic ) );
                b = lepton->spParent( il ) *
                    ( ( pslash * PR ) * charmquark.spParent( ic ) );
            }

            //EvtGenReport(EVTGEN_INFO,"EvtGen") <<"pslash*PR:"<<pslash*PR<<endl;
            //EvtGenReport(EVTGEN_INFO,"EvtGen") <<"sp charm:"<<charmquark.spParent(ic)<<endl;
            //EvtGenReport(EVTGEN_INFO,"EvtGen") <<"sp lepton:"<<lepton->spParent(il)<<endl;

            M[ic][il] = f1 * ( a1 * a + b1 * b ) + f2 * ( a2 * a + b2 * b );

            //EvtGenReport(EVTGEN_INFO,"EvtGen") << "Contr1:"<<a1<<" "<<a<<" "<<b1<<" "<<b<<endl;
            //EvtGenReport(EVTGEN_INFO,"EvtGen") << "Contr2:"<<a2<<" "<<a<<" "<<b2<<" "<<b<<endl;

            //EvtGenReport(EVTGEN_INFO,"EvtGen") <<"case1:"<<f1<<" "<<a1<<" "<<b1<<" "<<a<<" "<<b<<endl;
            //EvtGenReport(EVTGEN_INFO,"EvtGen") <<"case2:"<<f2<<" "<<a2<<" "<<b2<<" "<<a<<" "<<b<<endl;
        }
    }

    double prob = real( M[0][0] * conj( M[0][0] ) + M[1][0] * conj( M[1][0] ) +
                        M[0][1] * conj( M[0][1] ) + M[1][1] * conj( M[1][1] ) );

    //EvtGenReport(EVTGEN_INFO,"EvtGen") <<"prob:"<<prob<<endl;

    setProb( prob );

    return;
}
