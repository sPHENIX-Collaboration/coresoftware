////////////////////////////////////////////////////////////////////////////////
//
//  eASTProtonPhysics.hh  
//  Proton hadronic physics constructor for eASTPhysicsList 
//                                                     
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//                                                    
////////////////////////////////////////////////////////////////////////////////

#ifndef eASTProtonPhysics_h
#define eASTProtonPhysics_h 1

#include "G4VPhysicsConstructor.hh"

class G4TheoFSGenerator;
class G4FTFModel;
class G4ExcitedStringDecay;
class G4LundStringFragmentation;
class G4GeneratorPrecompoundInterface;


class eASTProtonPhysics: public G4VPhysicsConstructor
{
  public:
    eASTProtonPhysics();
    ~eASTProtonPhysics();

    virtual void ConstructProcess() override;
    virtual void ConstructParticle() override;

  private: 
    G4TheoFSGenerator* ftfp = nullptr;
    G4FTFModel* stringModel = nullptr;
    G4ExcitedStringDecay* stringDecay = nullptr;
    G4LundStringFragmentation* fragModel = nullptr;
    G4GeneratorPrecompoundInterface* preCompoundModel = nullptr;
 
};

#endif
