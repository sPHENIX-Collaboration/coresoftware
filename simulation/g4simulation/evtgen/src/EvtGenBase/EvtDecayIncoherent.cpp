
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

#include "EvtGenBase/EvtDecayIncoherent.hh"

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtReport.hh"

void EvtDecayIncoherent::makeDecay( EvtParticle* p, bool recursive )
{
    //initialize this the hard way..
    //Lange June 26, 2000
    for ( size_t i = 0; i < static_cast<unsigned int>( MAX_DAUG ); i++ ) {
        spinDensitySet[i] = 0;
    }

    _daugsDecayedByParentModel = false;

    decay( p );
    p->setDecayProb( 1.0 );

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
            //if (p->getDaug(i)->getNDaug()==0){
            //only do this if the user has not already set the
            //spin density matrix herself.
            //Lange June 26, 2000
            if ( isDaughterSpinDensitySet( i ) == 0 ) {
                p->getDaug( i )->setSpinDensityForward( rho );
            } else {
                //EvtGenReport(EVTGEN_INFO,"EvtGen") << "spinDensitymatrix already set!!!\n";
                EvtSpinDensity temp = p->getDaug( i )->getSpinDensityForward();
                //	EvtGenReport(EVTGEN_INFO,"EvtGen") <<temp<<endl;
            }
            //Now decay the daughter.  Really!
            p->getDaug( i )->decay();
        }
    }
}
