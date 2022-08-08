
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

#include "EvtGenBase/EvtDecayProb.hh"

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;

void EvtDecayProb::makeDecay( EvtParticle* p, bool recursive )
{
    int ntimes = 10000;

    double dummy;

    do {
        _weight = 1.0;
        _daugsDecayedByParentModel = false;

        decay( p );

        ntimes--;

        _prob = _prob / _weight;

        dummy = getProbMax( _prob ) * EvtRandom::Flat();
        p->setDecayProb( _prob / getProbMax( _prob ) );

    } while ( ntimes && ( _prob < dummy ) );

    if ( ntimes == 0 ) {
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
            << "Tried accept/reject:10000"
            << " times, and rejected all the times!" << endl;
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
            << "Is therefore accepting the last event!" << endl;
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
            << "Decay of particle:" << EvtPDL::name( p->getId() ).c_str()
            << "(channel:" << p->getChannel() << ") with mass " << p->mass()
            << endl;

        for ( size_t ii = 0; ii < p->getNDaug(); ii++ ) {
            EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
                << "Daughter " << ii << ":"
                << EvtPDL::name( p->getDaug( ii )->getId() ).c_str()
                << " with mass " << p->getDaug( ii )->mass() << endl;
        }
    }

    EvtSpinDensity rho;
    rho.setDiag( p->getSpinStates() );
    p->setSpinDensityBackward( rho );
    if ( getPHOTOS() || EvtRadCorr::alwaysRadCorr() ) {
        EvtRadCorr::doRadCorr( p );
    }

    if ( !recursive )
        return;

    //Now decay the daughters.
    if ( !daugsDecayedByParentModel() ) {
        for ( size_t i = 0; i < p->getNDaug(); i++ ) {
            //Need to set the spin density of the daughters to be
            //diagonal.
            rho.setDiag( p->getDaug( i )->getSpinStates() );
            p->getDaug( i )->setSpinDensityForward( rho );

            //Now decay the daughter.  Really!
            p->getDaug( i )->decay();
        }
    }
}
