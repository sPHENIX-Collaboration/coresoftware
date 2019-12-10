// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        ProtonPhysics.hh                                             //
//  Description: Proton hadronic physics constructor for EICPhysicsList       //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        22 June 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ProtonPhysics_h
#define ProtonPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;


class ProtonPhysics: public G4VPhysicsConstructor
{
  public:
    ProtonPhysics();
    ~ProtonPhysics();

    virtual void ConstructProcess() override;
    virtual void ConstructParticle() override;

  private: 
    G4TheoFSGenerator* ftfp;
    G4FTFModel* stringModel;
    G4ExcitedStringDecay* stringDecay;
    G4LundStringFragmentation* fragModel;
    G4GeneratorPrecompoundInterface* preCompoundModel;
 
};

#endif
