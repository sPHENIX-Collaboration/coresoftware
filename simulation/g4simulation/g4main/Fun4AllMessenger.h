#ifndef G4MAIN_FUN4ALLMESSENGER_H
#define G4MAIN_FUN4ALLMESSENGER_H

#include <Geant4/G4UImessenger.hh>

#include <Geant4/G4String.hh>       // for G4String

class Fun4AllServer;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIdirectory;

class Fun4AllMessenger: public G4UImessenger
{
public:
  Fun4AllMessenger(Fun4AllServer *ffa);
  ~Fun4AllMessenger() override;

void SetNewValue(G4UIcommand*, G4String) override;

private:

  G4UIdirectory*             m_Fun4AllDir;
  G4UIcmdWithAnInteger*      m_RunCmd;
  Fun4AllServer *se;
};

#endif
