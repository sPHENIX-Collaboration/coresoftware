// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        GammaLeptoNuclearPhysics.hh                                  //
//  Description: Gamma-nuclear, electro-nuclear and muon-nuclear physics      //
//                 constructor for EICPhysicsList                             //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        20 July 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef GammaLeptoNuclearPhysics_h
#define GammaLeptoNuclearPhysics_h 1

#include <Geant4/G4VPhysicsConstructor.hh>
#include <Geant4/G4GammaParticipants.hh>
#include <Geant4/G4QGSModel.hh>

class G4TheoFSGenerator;
class G4ExcitedStringDecay;
class G4QGSMFragmentation;
class G4GeneratorPrecompoundInterface;


class GammaLeptoNuclearPhysics: public G4VPhysicsConstructor
{
  public:
    GammaLeptoNuclearPhysics();
    ~GammaLeptoNuclearPhysics();

    virtual void ConstructProcess() override;
    virtual void ConstructParticle() override;

  private: 
    G4TheoFSGenerator* qgsp;
    G4QGSModel<G4GammaParticipants>* stringModel;
    G4ExcitedStringDecay* stringDecay;
    G4QGSMFragmentation* fragModel;
    G4GeneratorPrecompoundInterface* preCompoundModel;
 
};

#endif
