// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        IonPhysics.hh                                                //
//  Description: Ion hadronic physics constructor for EICPhysicsList          //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        6 July 2018                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef IonPhysics_h
#define IonPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;
class G4VComponentCrossSection;
class G4ComponentGGNuclNuclXsc;


class IonPhysics: public G4VPhysicsConstructor
{
  public:
    IonPhysics();
    ~IonPhysics();

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

    G4VComponentCrossSection* theGGNuclNuclXS; 
    G4ComponentGGNuclNuclXsc* ionGGXS;
};

#endif
