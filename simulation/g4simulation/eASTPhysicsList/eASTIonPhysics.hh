////////////////////////////////////////////////////////////////////////////////
// 
//  eASTIonPhysics.hh  
//  Ion hadronic physics constructor for eASTPhysicsList 
//                                                      
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//                                                     
////////////////////////////////////////////////////////////////////////////////

#ifndef eASTIonPhysics_h
#define eASTIonPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;
class G4VComponentCrossSection;
class G4ComponentGGNuclNuclXsc;


class eASTIonPhysics: public G4VPhysicsConstructor
{
  public:
    eASTIonPhysics();
    ~eASTIonPhysics();

    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    virtual void TerminateWorker() override;

  private: 
    G4TheoFSGenerator* ftfp = nullptr;
    G4FTFModel* stringModel = nullptr;
    G4ExcitedStringDecay* stringDecay = nullptr;
    G4LundStringFragmentation* fragModel = nullptr;
    G4GeneratorPrecompoundInterface* preCompoundModel = nullptr;

    G4VComponentCrossSection* theGGNuclNuclXS = nullptr; 
    G4ComponentGGNuclNuclXsc* ionGGXS = nullptr;
};

#endif
