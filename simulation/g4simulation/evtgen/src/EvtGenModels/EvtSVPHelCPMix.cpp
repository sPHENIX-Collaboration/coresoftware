
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

#include "EvtGenModels/EvtSVPHelCPMix.hh"

#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtSVPHelAmp.hh"

#include <iostream>
#include <stdlib.h>

std::string EvtSVPHelCPMix::getName()
{
    return "SVPHELCPMIX";
}

EvtDecayBase* EvtSVPHelCPMix::clone()
{
    return new EvtSVPHelCPMix;
}

void EvtSVPHelCPMix::init()
{
    // check that there are 5 arguments
    checkNArg( 5 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    checkSpinDaughter( 1, EvtSpinType::PHOTON );
}

void EvtSVPHelCPMix::initProbMax()
{
    setProbMax( 2.0 * ( getArg( 0 ) * getArg( 0 ) + getArg( 2 ) * getArg( 2 ) ) );
}

void EvtSVPHelCPMix::decay( EvtParticle* p )
{
    static EvtId BS0 = EvtPDL::getId( "B_s0" );
    //static EvtId BSB = EvtPDL::getId("anti-B_s0");

    //Flavour tagging of the initial state. Note that flavour mixing has already been applied out of this model
    //Initial_state == 0 (Bs at the initial state) and Initial_state == 1 (Anti-Bs in the initial state)
    int Initial_state( -1 );
    if ( EvtCPUtil::getInstance()->isBsMixed(
             p ) ) {    //The decaying particle has suffered flavour mixing, thus the initial state is its antiparticle
        if ( p->getId() == BS0 ) {
            Initial_state = 1;
        } else {
            Initial_state = 0;
        }
    } else {    //The decaying particle has NOT suffered flavour mixing, thus the initial state is itself
        if ( p->getId() == BS0 ) {
            Initial_state = 0;
        } else {
            Initial_state = 1;
        }
    }

    static EvtId BSH = EvtPDL::getId( "B_s0H" );
    static double ctauH = EvtPDL::getctau( BSH );
    static double gammaH = 1.0 / ctauH;

    static double deltaGamma = EvtCPUtil::getInstance()->getDeltaGamma( BS0 );

    //Here we're gonna generate and set the "envelope" lifetime, so we take the longest living component (for positive deltaGamma: tauH)
    //t is initialized following a e^(gammaH*t) lifetime distribution. When computing the amplitudes a factor e^(gammaH*t/2) should be substracted.
    double t = -log( EvtRandom::Flat() ) *
               ( 1.0 / gammaH );    //This overrules the lifetimes made by the program performing the mixing (CPUtil)
    if ( EvtCPUtil::getInstance()->isBsMixed( p ) ) {
        p->getParent()->setLifetime( t );
    } else {
        p->setLifetime( t );
    }

    static double deltaMs = EvtCPUtil::getInstance()->getDeltaM( BS0 );
    double mt = exp( -std::max( 0.0, deltaGamma ) * t / ( 2.0 * EvtConst::c ) );
    double pt = exp( +std::min( 0.0, deltaGamma ) * t / ( 2.0 * EvtConst::c ) );

    //Using the same sign convention as in J.P. Silva, hep-ph/0410351 (2004)
    EvtComplex qp = EvtComplex( cos( -2.0 * getArg( 4 ) ),
                                sin( -2.0 * getArg( 4 ) ) );    // q/p=e^(-2*beta_s)
    EvtComplex gplus =
        ( mt * EvtComplex( cos( deltaMs * t / ( 2.0 * EvtConst::c ) ),
                           sin( deltaMs * t / ( 2.0 * EvtConst::c ) ) ) +
          pt * EvtComplex( cos( deltaMs * t / ( 2.0 * EvtConst::c ) ),
                           sin( -deltaMs * t / ( 2.0 * EvtConst::c ) ) ) ) /
        2.0;
    EvtComplex gminus =
        ( +mt * EvtComplex( cos( deltaMs * t / ( 2.0 * EvtConst::c ) ),
                            sin( deltaMs * t / ( 2.0 * EvtConst::c ) ) ) -
          pt * EvtComplex( cos( deltaMs * t / ( 2.0 * EvtConst::c ) ),
                           sin( -deltaMs * t / ( 2.0 * EvtConst::c ) ) ) ) /
        2.0;

    //These should be filled with the helicity amplitudes at t=0
    EvtComplex arg_hm, arg_hp;
    arg_hp = EvtComplex( getArg( 0 ) * cos( getArg( 1 ) ),
                         getArg( 0 ) * sin( getArg( 1 ) ) );
    arg_hm = EvtComplex( getArg( 2 ) * cos( getArg( 3 ) ),
                         getArg( 2 ) * sin( getArg( 3 ) ) );

    //Time-dependent amplitudes H+(t) and H-(t) are computed for a Bs and Anti-Bs in the initial state
    EvtComplex hp, hm;
    if ( Initial_state == 0 ) {    //These are the equations for Bs

        hp = arg_hp * gplus + qp * conj( arg_hm ) * gminus;
        hm = arg_hm * gplus + qp * conj( arg_hp ) * gminus;

    } else if ( Initial_state == 1 ) {    //The equations for Anti-Bs

        hp = conj( arg_hm ) * gplus + ( 1.0 / qp ) * arg_hp * gminus;
        hm = conj( arg_hp ) * gplus + ( 1.0 / qp ) * arg_hm * gminus;

    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Initial state was not BSB or BS0!" << std::endl;
        ::abort();
    }

    //Compute the decay amplitudes from the time-dependent helicity amplitudes
    EvtSVPHelAmp::SVPHel( p, _amp2, getDaug( 0 ), getDaug( 1 ), hp, hm );

    return;
}
