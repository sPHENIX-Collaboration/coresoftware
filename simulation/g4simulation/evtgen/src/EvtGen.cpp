
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

#include "EvtGen/EvtGen.hh"

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtStatus.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include "EvtGenModels/EvtModelReg.hh"
#include "EvtGenModels/EvtNoRadCorr.hh"

#include <cstdlib>
#include <fstream>
#include <string>

using std::endl;
using std::fstream;
using std::ifstream;

EvtGen::~EvtGen()
{
    //This is a bit ugly, should not do anything
    //in a destructor. This will fail if EvtGen is made a static
    //because then this destructor might be called _after_
    //the destruction of objects that it depends on, e.g., EvtPDL.

    if ( getenv( "EVTINFO" ) ) {
        EvtDecayTable::getInstance()->printSummary();
    }
}

EvtGen::EvtGen( const std::string& decayName, const std::string& pdtTableName,
                EvtRandomEngine* randomEngine, EvtAbsRadCorr* isrEngine,
                const std::list<EvtDecayBase*>* extraModels, int mixingType,
                bool useXml )
{
    std::ifstream pdtIn( pdtTableName );
    if ( !pdtIn ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Could not open:" << pdtTableName << "EvtPDL" << endl;
        return;
    }
    initialize( decayName, pdtIn, randomEngine, isrEngine, extraModels,
                mixingType, useXml );
    pdtIn.close();
}

EvtGen::EvtGen( const std::string& decayName, std::istream& pdtTableData,
                EvtRandomEngine* randomEngine, EvtAbsRadCorr* isrEngine,
                const std::list<EvtDecayBase*>* extraModels, int mixingType,
                bool useXml )
{
    initialize( decayName, pdtTableData, randomEngine, isrEngine, extraModels,
                mixingType, useXml );
}

void EvtGen::initialize( const std::string& decayName, std::istream& pdtTable,
                         EvtRandomEngine* randomEngine, EvtAbsRadCorr* isrEngine,
                         const std::list<EvtDecayBase*>* extraModels,
                         int mixingType, bool useXml )
{
    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Initializing EvtGen" << endl;

    if ( randomEngine == 0 ) {
        static EvtSimpleRandomEngine defaultRandomEngine;
        EvtRandom::setRandomEngine( &defaultRandomEngine );
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "No random engine given in "
            << "EvtGen::EvtGen constructor, "
            << "will use default EvtSimpleRandomEngine." << endl;
    } else {
        EvtRandom::setRandomEngine( randomEngine );
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Storing known decay models" << endl;
    EvtModelReg dummy( extraModels );

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Main decay file name  :" << decayName << endl;

    _pdl.readPDT( pdtTable );

    if ( useXml ) {
        EvtDecayTable::getInstance()->readXMLDecayFile( decayName, false );
    } else {
        EvtDecayTable::getInstance()->readDecayFile( decayName, false );
    }

    _mixingType = mixingType;
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Mixing type integer set to " << _mixingType << endl;
    EvtCPUtil::getInstance()->setMixingType( _mixingType );

    // Set the radiative correction engine

    if ( isrEngine != 0 ) {
        EvtRadCorr::setRadCorrEngine( isrEngine );

    } else {
        // Owing to the pure abstract interface, we still need to define a concrete
        // implementation of a radiative correction engine. Use one which does nothing.
        EvtAbsRadCorr* noRadCorr = new EvtNoRadCorr();
        EvtRadCorr::setRadCorrEngine( noRadCorr );
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Done initializing EvtGen" << endl;
}

void EvtGen::readUDecay( const std::string& uDecayName, bool useXml )
{
    ifstream indec;

    if ( uDecayName.size() == 0 ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Is not reading a user decay file!" << endl;
    } else {
        indec.open( uDecayName );
        if ( indec ) {
            if ( useXml ) {
                EvtDecayTable::getInstance()->readXMLDecayFile( uDecayName, true );
            } else {
                EvtDecayTable::getInstance()->readDecayFile( uDecayName, true );
            }
        } else {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "Can not find UDECAY file '" << uDecayName << "'.  Stopping"
                << endl;
            ::abort();
        }
    }
}

EvtHepMCEvent* EvtGen::generateDecay( int PDGId, EvtVector4R refFrameP4,
                                      EvtVector4R translation,
                                      EvtSpinDensity* spinDensity )
{
    EvtParticle* theParticle( 0 );

    if ( spinDensity == 0 ) {
        theParticle = EvtParticleFactory::particleFactory(
            EvtPDL::evtIdFromStdHep( PDGId ), refFrameP4 );
    } else {
        theParticle = EvtParticleFactory::particleFactory(
            EvtPDL::evtIdFromStdHep( PDGId ), refFrameP4, *spinDensity );
    }

    generateDecay( theParticle );
    EvtHepMCEvent* hepMCEvent = new EvtHepMCEvent();
    hepMCEvent->constructEvent( theParticle, translation );

    theParticle->deleteTree();

    return hepMCEvent;
}

void EvtGen::generateDecay( EvtParticle* p )
{
    int times = 0;
    do {
        times += 1;
        EvtStatus::initRejectFlag();

        p->decay();
        //ok then finish.
        if ( EvtStatus::getRejectFlag() == 0 ) {
            times = 0;
        } else {
            for ( size_t ii = 0; ii < p->getNDaug(); ii++ ) {
                EvtParticle* temp = p->getDaug( ii );
                temp->deleteTree();
            }
            p->resetFirstOrNot();
            p->resetNDaug();
        }

        if ( times == 10000 ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Your event has been rejected 10000 times!" << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Will now abort." << endl;
            ::abort();
            times = 0;
        }
    } while ( times );
}
