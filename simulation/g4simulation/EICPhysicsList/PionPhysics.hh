// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        PionPhysics.hh                                               //
//  Description: Pion hadronic physics constructor for EICPhysicsList         //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        21 June 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef PionPhysics_h
#define PionPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;


class PionPhysics: public G4VPhysicsConstructor
{
  public:
    PionPhysics();
    ~PionPhysics();

    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;

  private: 
    G4TheoFSGenerator* ftfp;
    G4FTFModel* stringModel;
    G4ExcitedStringDecay* stringDecay;
    G4LundStringFragmentation* fragModel;
    G4GeneratorPrecompoundInterface* preCompoundModel;
 
};

#endif
