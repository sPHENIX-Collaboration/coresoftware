// sPHENIX prototype detector messager
// Dec 2013, hexc
// Jan 2, 2014:  Added more messenger commands

#include "PHG4HcalPrototypeDetectorMessenger.h"

#include "PHG4HcalPrototypeDetector.h"

#include <Geant4/G4ApplicationState.hh>         // for G4State_Idle, G4State...
#include <Geant4/G4UIdirectory.hh>
#include <Geant4/G4UIcmdWithAString.hh>
#include <Geant4/G4UIcmdWithADoubleAndUnit.hh>
#include <Geant4/G4UIcmdWithoutParameter.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PHG4HcalPrototypeDetectorMessenger::PHG4HcalPrototypeDetectorMessenger(PHG4HcalPrototypeDetector * Det)
:fPHG4HcalPrototypeDetector(Det)
{ 
  sPhnxDir = new G4UIdirectory("/sPhnx/");
  sPhnxDir->SetGuidance("UI commands for modifying the prototype detector properties.");

  detDir = new G4UIdirectory("/sPhnx/det/");
  detDir->SetGuidance("UI commands for changing the detector geometry.");

  fMaterCmd = new G4UIcmdWithAString("/sPhnx/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the world.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
       
  fUpdateCmd = new G4UIcmdWithoutParameter("/sPhnx/det/update",this);
  fUpdateCmd->SetGuidance("Update geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

  outerHcalDPhiCmd = new G4UIcmdWithADoubleAndUnit("/sPhnx/det/setOuterHcalAngularSpan",this);  
  outerHcalDPhiCmd->SetGuidance("Change the outer hcal angular span (default: 0.28 rad).");
  outerHcalDPhiCmd->SetGuidance("Need to run the detector update command after running this command!");
  outerHcalDPhiCmd->SetParameterName("hcal2DPhi",false);
  outerHcalDPhiCmd->SetUnitCategory("radian");
  outerHcalDPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  outerPlateTiltAngleCmd = new G4UIcmdWithADoubleAndUnit("/sPhnx/det/setOuterPlateTiltAngule",this);  
  outerPlateTiltAngleCmd->SetGuidance("Change the outer plate tilt angle (default: 0.12 rad).");
  outerPlateTiltAngleCmd->SetGuidance("Need to run the detector update command after running this command!");
  outerPlateTiltAngleCmd->SetParameterName("hcal2TiltAngle",false);
  outerPlateTiltAngleCmd->SetUnitCategory("radian");
  outerPlateTiltAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  innerHcalDPhiCmd = new G4UIcmdWithADoubleAndUnit("/sPhnx/det/setInnerHcalAngularSpan",this);  
  innerHcalDPhiCmd->SetGuidance("Change the inner hcal angular span (default: 0.28 rad).");
  innerHcalDPhiCmd->SetGuidance("Need to run the detector update command after running this command!");
  innerHcalDPhiCmd->SetParameterName("hcal1DPhi",false);
  innerHcalDPhiCmd->SetUnitCategory("radian");
  innerHcalDPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  innerPlateTiltAngleCmd = new G4UIcmdWithADoubleAndUnit("/sPhnx/det/setInnerPlateTiltAngule",this);  
  innerPlateTiltAngleCmd->SetGuidance("Change the inner plate tilt angle (default: 0.30 rad).");
  innerPlateTiltAngleCmd->SetGuidance("Need to run the detector update command after running this command!");
  innerPlateTiltAngleCmd->SetParameterName("hcal1TiltAngle",false);
  innerPlateTiltAngleCmd->SetUnitCategory("radian");
  innerPlateTiltAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PHG4HcalPrototypeDetectorMessenger::~PHG4HcalPrototypeDetectorMessenger()
{
  delete fMaterCmd;
  delete fUpdateCmd;
  delete outerHcalDPhiCmd;
  delete innerHcalDPhiCmd;
  delete outerPlateTiltAngleCmd;
  delete innerPlateTiltAngleCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PHG4HcalPrototypeDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fPHG4HcalPrototypeDetector->SetMaterial(newValue);}
   
  if( command == fUpdateCmd )
   { fPHG4HcalPrototypeDetector->UpdateGeometry();}

  if ( command == outerHcalDPhiCmd )
    { fPHG4HcalPrototypeDetector->SetOuterHcalDPhi(outerHcalDPhiCmd->GetNewDoubleValue(newValue));}

  if ( command == outerPlateTiltAngleCmd )
    { fPHG4HcalPrototypeDetector->SetOuterPlateTiltAngle(outerPlateTiltAngleCmd->GetNewDoubleValue(newValue));}

  if ( command == innerHcalDPhiCmd )
    { fPHG4HcalPrototypeDetector->SetInnerHcalDPhi(innerHcalDPhiCmd->GetNewDoubleValue(newValue));}

  if ( command == innerPlateTiltAngleCmd )
    { fPHG4HcalPrototypeDetector->SetInnerPlateTiltAngle(innerPlateTiltAngleCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
