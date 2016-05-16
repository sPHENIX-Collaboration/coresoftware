// header file for the 2nd HCal prototype detector construction messager
// Nov, 2015 - hexc

#ifndef PHG4HcalPrototype2DetectorMessenger_h
#define PHG4HcalPrototype2DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PHG4HcalPrototype2Detector;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PHG4HcalPrototype2DetectorMessenger: public G4UImessenger
{
public:
  
  PHG4HcalPrototype2DetectorMessenger(PHG4HcalPrototype2Detector* );
  ~PHG4HcalPrototype2DetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  PHG4HcalPrototype2Detector*      fPHG4HcalPrototype2Detector;
  
  G4UIdirectory*             sPhnxDir;
  G4UIdirectory*             detDir;

  G4UIcmdWithAString*        fMaterCmd;
  G4UIcmdWithoutParameter*   fUpdateCmd;
  G4UIcmdWithADoubleAndUnit  *outerHcalDPhiCmd, *outerPlateTiltAngleCmd;
  G4UIcmdWithADoubleAndUnit  *innerHcalDPhiCmd, *innerPlateTiltAngleCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

