
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

#ifndef EVTHEPMCEVENT_HH
#define EVTHEPMCEVENT_HH

#include "EvtGenBase/EvtVector4R.hh"

#ifdef EVTGEN_HEPMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/Units.h"
typedef HepMC3::GenParticlePtr GenParticlePtr;
typedef HepMC3::GenVertexPtr GenVertexPtr;
typedef HepMC3::GenEvent GenEvent;
typedef HepMC3::FourVector FourVector;
typedef HepMC3::Units Units;
inline GenParticlePtr newGenParticlePtr(
    const FourVector& mom = FourVector::ZERO_VECTOR(), int pid = 0,
    int status = 0 )
{
    std::cout << "USING HEPMC 3 BRO" << std::endl;
    return std::make_shared<HepMC3::GenParticle>( mom, pid, status );
}
inline GenVertexPtr newGenVertexPtr(
    const FourVector& pos = FourVector::ZERO_VECTOR() )
{
    std::cout << "USING HEPMC 3 BRO" << std::endl;	
    return std::make_shared<HepMC3::GenVertex>( pos );
}
#else
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/Units.h"
typedef HepMC::GenParticle* GenParticlePtr;
typedef HepMC::GenVertex* GenVertexPtr;
typedef HepMC::GenEvent GenEvent;
typedef HepMC::FourVector FourVector;
#define Units HepMC::Units
inline GenParticlePtr newGenParticlePtr(
    const FourVector& mom = FourVector( 0.0, 0.0, 0.0, 0.0 ), int pid = 0,
    int status = 0 )
{
  //  std::cout << "USING JUST HEPMC BRO" << std::endl;	
    return new HepMC::GenParticle( mom, pid, status );
}
inline GenVertexPtr newGenVertexPtr(
    const FourVector& pos = FourVector( 0.0, 0.0, 0.0, 0.0 ) )
{
   // std::cout << "USING JUST HEPMC BRO" << std::endl;		
    return new HepMC::GenVertex( pos );
}
#endif

class EvtParticle;

class EvtHepMCEvent {
  public:
    EvtHepMCEvent();
    virtual ~EvtHepMCEvent();

    // Select what frame a given GenParticle is in:
    // its own restframe, the lab frame (first mother), or its mother's frame
    enum HepMCFrame
    {
        RESTFRAME = 1,
        LAB = 2,
        MOTHER = 3
    };
    // Select the GenParticle status
    enum HepMCStatus
    {
        STABLE = 1,
        DECAYED = 2,
        HISTORY = 3
    };

    void constructEvent( EvtParticle* baseParticle );
    void constructEvent( EvtParticle* baseParticle, EvtVector4R& translation );

    GenEvent* getEvent() { return _theEvent; }

    // Methods used to create GenParticles and FourVectors of vertices.
    // Make these public so that other classes may call them if they use EvtHepMCEvent.

    // Create a GenParticle using info from the EvtParticle, specifying what frame
    // the 4-momentum is from.
    GenParticlePtr createGenParticle( EvtParticle* theParticle, int frameType );

    // Find out the decay vertex position for the given EvtParticle.
    FourVector getVertexCoord( EvtParticle* theParticle );

  protected:
  private:
    // Delete the event structure (called by destructor)
    void deleteEvent();

    // Add a vertex to the event. This is called by the constructEvent function
    // and is recursive, i.e. it loops through all possible daughter particles and
    // their descendents.
    void addVertex( EvtParticle* inEvtParticle, GenParticlePtr inGenParticle );

    GenEvent* _theEvent;
    EvtVector4R _translation;
};

#endif
