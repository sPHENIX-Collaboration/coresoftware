// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        AntiBaryonPhysics.hh                                         //
//  Description: Anti-baryon hadronic physics constructor for EICPhysicsList  //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        5 July 2018                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AntiBaryonPhysics_h
#define AntiBaryonPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;
class G4ComponentAntiNuclNuclearXS; 

class AntiBaryonPhysics: public G4VPhysicsConstructor
{
  public:
    AntiBaryonPhysics();
    ~AntiBaryonPhysics();

    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
#if G4VERSION_NUMBER >= 1004
  virtual void TerminateWorker() override {}
#endif

  private: 
    G4TheoFSGenerator* ftfp;
    G4FTFModel* stringModel;
    G4ExcitedStringDecay* stringDecay;
    G4LundStringFragmentation* fragModel;
    G4GeneratorPrecompoundInterface* preCompoundModel;

    G4ComponentAntiNuclNuclearXS* theAntiNucleonXS; 
};

#endif
