// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        KaonPhysics.hh                                               //
//  Description: Kaon hadronic physics constructor for EICPhysicsList         //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        3 July 2018                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef KaonPhysics_h
#define KaonPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;


class KaonPhysics: public G4VPhysicsConstructor
{
  public:
    KaonPhysics();
    ~KaonPhysics();

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
 
};

#endif
