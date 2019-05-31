// header file for the detector construction messager
// Jan 2, 2014 hexc

#ifndef G4DETECTORS_PHG4HCALPROTOTYPEDETECTORMESSENGER_H
#define G4DETECTORS_PHG4HCALPROTOTYPEDETECTORMESSENGER_H

#include <Geant4/globals.hh>
#include <Geant4/G4UImessenger.hh>

class PHG4HcalPrototypeDetector;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PHG4HcalPrototypeDetectorMessenger: public G4UImessenger
{
public:
  
  PHG4HcalPrototypeDetectorMessenger(PHG4HcalPrototypeDetector* );
  ~PHG4HcalPrototypeDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  PHG4HcalPrototypeDetector*      fPHG4HcalPrototypeDetector;
  
  G4UIdirectory*             sPhnxDir;
  G4UIdirectory*             detDir;

  G4UIcmdWithAString*        fMaterCmd;
  G4UIcmdWithoutParameter*   fUpdateCmd;
  G4UIcmdWithADoubleAndUnit  *outerHcalDPhiCmd, *outerPlateTiltAngleCmd;
  G4UIcmdWithADoubleAndUnit  *innerHcalDPhiCmd, *innerPlateTiltAngleCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

