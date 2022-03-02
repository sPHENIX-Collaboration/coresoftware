// ********************************************************************
//
// eASTPhysicsListMessenger.hh
//   Header file of a messenger class that handles physics options.
//
// History
//   May 8th, 2021 : first implementation
//
// ********************************************************************

#ifndef eASTPhysicsListMessenger_H
#define eASTPhysicsListMessenger_H 1

#include "G4UImessenger.hh"
#include "globals.hh"

class eASTPhysicsList;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class eASTPhysicsListMessenger: public G4UImessenger
{
  public:
    eASTPhysicsListMessenger(eASTPhysicsList*);
    virtual ~eASTPhysicsListMessenger();
    virtual void SetNewValue(G4UIcommand*,G4String);
    virtual G4String GetCurrentValue(G4UIcommand*);

  private:
    eASTPhysicsList* pPL;
    G4UIdirectory* physDir;
    //G4UIcmdWithoutParameter*   addHPCmd;
    G4UIcmdWithoutParameter*   addRDMCmd;
    G4UIcmdWithoutParameter*   addOpticalCmd;
    G4UIcmdWithAString* addStepLimitCmd;

    G4UIdirectory* physLimitDir;
    G4UIcmdWithADoubleAndUnit* setStepLimitCmd;
    G4UIcommand*        setRegionStepLimitCmd;

    G4UIdirectory* physCutDir;
    G4UIcmdWithADoubleAndUnit* setCutCmd;
    G4UIcommand*        setCutParticleCmd;
    G4UIcommand*        setCutRegionCmd;
    G4UIcommand*        setCutRegionParticleCmd;

};

#endif

