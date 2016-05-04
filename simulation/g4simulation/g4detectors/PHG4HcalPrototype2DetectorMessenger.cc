// sPHENIX 2nd prototype detector messager
// Nov 2015, hexc


#include "PHG4HcalPrototype2DetectorMessenger.h"

#include "PHG4HcalPrototype2Detector.h"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PHG4HcalPrototype2DetectorMessenger::PHG4HcalPrototype2DetectorMessenger(PHG4HcalPrototype2Detector * Det)
:fPHG4HcalPrototype2Detector(Det)
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

PHG4HcalPrototype2DetectorMessenger::~PHG4HcalPrototype2DetectorMessenger()
{
  delete fMaterCmd;
  delete fUpdateCmd;
  delete outerHcalDPhiCmd;
  delete innerHcalDPhiCmd;
  delete outerPlateTiltAngleCmd;
  delete innerPlateTiltAngleCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PHG4HcalPrototype2DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fPHG4HcalPrototype2Detector->SetMaterial(newValue);}
   
  if( command == fUpdateCmd )
   { fPHG4HcalPrototype2Detector->UpdateGeometry();}

  if ( command == outerHcalDPhiCmd )
    { fPHG4HcalPrototype2Detector->SetOuterHcalDPhi(outerHcalDPhiCmd->GetNewDoubleValue(newValue));}

  if ( command == outerPlateTiltAngleCmd )
    { fPHG4HcalPrototype2Detector->SetOuterPlateTiltAngle(outerPlateTiltAngleCmd->GetNewDoubleValue(newValue));}

  if ( command == innerHcalDPhiCmd )
    { fPHG4HcalPrototype2Detector->SetInnerHcalDPhi(innerHcalDPhiCmd->GetNewDoubleValue(newValue));}

  if ( command == innerPlateTiltAngleCmd )
    { fPHG4HcalPrototype2Detector->SetInnerPlateTiltAngle(innerPlateTiltAngleCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
