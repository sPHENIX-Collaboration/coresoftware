
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

#include "EvtGenBase/EvtHepMCEvent.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"

EvtHepMCEvent::EvtHepMCEvent() :
    _theEvent( 0 ), _translation( 0.0, 0.0, 0.0, 0.0 )
{
}

EvtHepMCEvent::~EvtHepMCEvent()
{
    this->deleteEvent();
}

void EvtHepMCEvent::deleteEvent()
{
    if ( _theEvent != 0 ) {
        _theEvent->clear();
        delete _theEvent;
        _theEvent = 0;
    }
}

void EvtHepMCEvent::constructEvent( EvtParticle* baseParticle )
{
    EvtVector4R origin( 0.0, 0.0, 0.0, 0.0 );
    this->constructEvent( baseParticle, origin );
	//std::cout << "Pass 1 Here inside EvtHepMCEvent ORB ZZ ROCK" << std::endl;
}

void EvtHepMCEvent::constructEvent( EvtParticle* baseParticle,
                                    EvtVector4R& translation )
{
    // This class does not take ownership of the base particle pointer.
    // Rather, it uses the base particle to construct the event.
//	std::cout << "Pass 2 Here inside EvtHepMCEvent ORB ZZ ROCK" << std::endl;

    this->deleteEvent();
    if ( baseParticle == 0 ) {
	//	std::cout << "No FUCKIN BASE Particles BRO" << std::endl;

        return;
    }
//	std::cout << "Pass 3 Here inside EvtHepMCEvent ORB ZZ ROCK" << std::endl;

    _theEvent = new GenEvent( Units::GEV, Units::MM );
    _translation = translation;

    // Use the recursive function addVertex to add a vertex with incoming/outgoing
    // particles. Adds a new vertex for any EvtParticles with decay daughters.
    // All particles are in the rest frame of the base particle ("lab frame").
	//std::cout << "Pass 4 Here inside EvtHepMCEvent ORB ZZ ROCK" << std::endl;

    GenParticlePtr hepMCGenParticle =
        this->createGenParticle( baseParticle, EvtHepMCEvent::LAB );

    this->addVertex( baseParticle, hepMCGenParticle );
	//std::cout << "Pass 5 Here inside EvtHepMCEvent ORB ZZ ROCK" << std::endl;
	
}

GenParticlePtr EvtHepMCEvent::createGenParticle( EvtParticle* theParticle,
                                                 int frameType )
{
//	std::cout << "Pass 1 Here inside createGenParticle ORB ZZ ROCK" << std::endl;


    // Create an HepMC GenParticle, with the 4-momenta in the frame given by the frameType integer
    GenParticlePtr genParticle{nullptr};
//	std::cout << "Pass 2 Here inside createGenParticle ORB ZZ ROCK" << std::endl;

    if ( theParticle != 0 ) {
        // Set the particle status integer to either stable or decayed
        int status( EvtHepMCEvent::STABLE );
        int nDaug = theParticle->getNDaug();
        if ( nDaug > 0 ) {
            status = EvtHepMCEvent::DECAYED;
        }

        // Get the 4-momentum (E, px, py, pz) for the EvtParticle.
        EvtVector4R p4( 0.0, 0.0, 0.0, 0.0 );
//		std::cout << "Pass 2.5 Here inside createGenParticle ORB ZZ ROCK" << std::endl;

        if ( frameType == EvtHepMCEvent::RESTFRAME ) {
            p4 = theParticle->getP4Restframe();
        } else if ( frameType == EvtHepMCEvent::LAB ) {
            p4 = theParticle->getP4Lab();
        } else {
            p4 = theParticle->getP4();
        }
//		std::cout << "Pass 2.7 Here inside createGenParticle ORB ZZ ROCK" << std::endl;

        // Convert this to the HepMC 4-momentum
        double E = p4.get( 0 );
        double px = p4.get( 1 );
        double py = p4.get( 2 );
        double pz = p4.get( 3 );

        FourVector hepMC_p4( px, py, pz, E );

        // Get the particle PDG integer id
        int PDGInt = EvtPDL::getStdHep( theParticle->getId() );

        genParticle = newGenParticlePtr( hepMC_p4, PDGInt, status );
//	std::cout << "Pass 3 Here inside createGenParticle ORB ZZ ROCK" << std::endl;
		
    }
	//std::cout << "Pass 4 Here inside createGenParticle ORB ZZ ROCK" << std::endl;

    return genParticle;
}

void EvtHepMCEvent::addVertex( EvtParticle* inEvtParticle,
                               GenParticlePtr inGenParticle )
{
    // This is a recursive function that adds GenVertices to the GenEvent for
    // the incoming EvtParticle and its daughters. We use two separate
    // pointers for the EvtParticle and GenParticle information: the former
    // to obtain the PDGId, 4-momenta, daughter and vertex positions, the latter to
    // set the incoming particle to the vertex. Note that the outgoing particle for
    // one vertex might be the incoming particle for another vertex - this needs to
    // be the same GenParticle pointer, hence the reason for using it as a 2nd argument
    // in this function.
	//std::cout << "Now Add Vertex Pass 1" << std::endl;

    if ( _theEvent == 0 || inEvtParticle == 0 || inGenParticle == 0 ) {
        return;
	//	std::cout << "Now Add Vertex Pass 2" << std::endl;
    }

    // Create the decay vertex
    FourVector vtxCoord = this->getVertexCoord( inEvtParticle );
    GenVertexPtr theVertex = newGenVertexPtr( vtxCoord );

//	std::cout << "Now Add Vertex Pass 3" << std::endl;


    // Add the vertex to the event
    _theEvent->add_vertex( theVertex );

    // Set the incoming particle
    theVertex->add_particle_in( inGenParticle );

    // Set the outgoing particles (decay products)
    int nDaug = inEvtParticle->getNDaug();
    int iDaug( 0 );
//	std::cout << "Now Add Vertex Pass 4 and nDaug = " << nDaug  << std::endl;
	
    // Loop over the daughters
    for ( iDaug = 0; iDaug < nDaug; iDaug++ ) {
        EvtParticle* evtDaughter = inEvtParticle->getDaug( iDaug );
        GenParticlePtr genDaughter =
            this->createGenParticle( evtDaughter, EvtHepMCEvent::LAB );

        if ( genDaughter != 0 ) {
            // Add a new GenParticle (outgoing) particle daughter to the vertex
            theVertex->add_particle_out( genDaughter );

            // Find out if the daughter also has decay products.
            // If so, recursively run this function again.
            int nDaugProducts = evtDaughter->getNDaug();

            if ( nDaugProducts > 0 ) {
                // Recursively process daughter particles and add their vertices to the event
                this->addVertex( evtDaughter, genDaughter );

            }    // Have daughter products

        }    // hepMCDaughter != 0

    }    // Loop over daughters

//	std::cout << "Now Add Vertex Pass 5 "  << std::endl;
	
}

FourVector EvtHepMCEvent::getVertexCoord( EvtParticle* theParticle )
{
    FourVector vertexCoord( 0.0, 0.0, 0.0, 0.0 );

    if ( theParticle != 0 && theParticle->getNDaug() != 0 ) {
        // Get the position (t,x,y,z) of the EvtParticle, offset by the translation vector.
        // This position will be the point where the particle decays. So we ask
        // the position of the (1st) daughter particle.
        EvtParticle* daugParticle = theParticle->getDaug( 0 );

        if ( daugParticle != 0 ) {
            EvtVector4R vtxPosition = daugParticle->get4Pos() + _translation;

            // Create the HepMC 4 vector of the position (x,y,z,t)
            vertexCoord.setX( vtxPosition.get( 1 ) );
            vertexCoord.setY( vtxPosition.get( 2 ) );
            vertexCoord.setZ( vtxPosition.get( 3 ) );
            vertexCoord.setT( vtxPosition.get( 0 ) );
        }
    }

    return vertexCoord;
}
