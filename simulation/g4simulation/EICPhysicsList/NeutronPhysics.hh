// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        NeutronPhysics.hh                                            //
//  Description: Neutron hadronic physics constructor for EICPhysicsList      //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        3 July 2018                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef NeutronPhysics_h
#define NeutronPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;


class NeutronPhysics: public G4VPhysicsConstructor
{
  public:
    NeutronPhysics();
    ~NeutronPhysics();

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
