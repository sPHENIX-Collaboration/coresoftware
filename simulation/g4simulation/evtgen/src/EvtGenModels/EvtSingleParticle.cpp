
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

#include "EvtGenModels/EvtSingleParticle.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtSingleParticle::getName()
{
    return "SINGLE";
}

EvtDecayBase* EvtSingleParticle::clone()
{
    return new EvtSingleParticle();
}

void EvtSingleParticle::init()
{
    //turn off checks for charge conservation
    disableCheckQ();

    if ( ( getNArg() == 6 ) || ( getNArg() == 4 ) || ( getNArg() == 2 ) ) {
        if ( getNArg() == 6 ) {
            //copy the arguments into eaiser to remember names!

            pmin = getArg( 0 );
            pmax = getArg( 1 );

            cthetamin = getArg( 2 );
            cthetamax = getArg( 3 );

            phimin = getArg( 4 );
            phimax = getArg( 5 );
        }

        if ( getNArg() == 4 ) {
            //copy the arguments into eaiser to remember names!

            pmin = getArg( 0 );
            pmax = getArg( 1 );

            cthetamin = getArg( 2 );
            cthetamax = getArg( 3 );

            phimin = 0.0;
            phimax = EvtConst::twoPi;
        }

        if ( getNArg() == 2 ) {
            //copy the arguments into eaiser to remember names!

            pmin = getArg( 0 );
            pmax = getArg( 1 );

            cthetamin = -1.0;
            cthetamax = 1.0;

            phimin = 0.0;
            phimax = EvtConst::twoPi;
        }

    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtSingleParticle generator expected "
            << " 6, 4, or 2 arguments but found:" << getNArg() << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "The single particle generator has been configured:" << endl;
    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << pmax << " > p > " << pmin << endl;
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << cthetamax << " > costheta > " << cthetamin << endl;
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << phimax << " > phi > " << phimin << endl;
}

void EvtSingleParticle::decay( EvtParticle* p )
{
    EvtParticle* d;
    EvtVector4R p4;

    double mass = EvtPDL::getMass( getDaug( 0 ) );

    p->makeDaughters( getNDaug(), getDaugs() );
    d = p->getDaug( 0 );

    //generate flat distribution in p
    //we are now in the parents restframe! This means the
    //restframe of the e+e- collison.
    double pcm = EvtRandom::Flat( pmin, pmax );
    //generate flat distribution in phi.
    double phi = EvtRandom::Flat( phimin, phimax );

    double cthetalab;

    do {
        //generate flat distribution in costheta
        double ctheta = EvtRandom::Flat( cthetamin, cthetamax );
        double stheta = sqrt( 1.0 - ctheta * ctheta );
        p4.set( sqrt( mass * mass + pcm * pcm ), pcm * cos( phi ) * stheta,
                pcm * sin( phi ) * stheta, pcm * ctheta );

        d->init( getDaug( 0 ), p4 );

        //get 4 vector in the lab frame!
        EvtVector4R p4lab = d->getP4Lab();
        cthetalab = p4lab.get( 3 ) / p4lab.d3mag();
    } while ( cthetalab > cthetamax || cthetalab < cthetamin );

    return;
}
