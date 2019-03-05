// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        EICPhysicsList.hh                                            //
//  Description: Geant4 physics list for Electron Ion Collider detectors      //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        21 June 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef EICPhysicsList_h
#define EICPhysicsList_h 1

#include <Geant4/G4VModularPhysicsList.hh>
#include <Geant4/globals.hh>


class EICPhysicsList: public G4VModularPhysicsList
{
  public:
    EICPhysicsList();
    ~EICPhysicsList();

    virtual void ConstructParticle();
    virtual void SetCuts();

};

#endif
