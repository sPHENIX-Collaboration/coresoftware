#include "Fun4AllMessenger.h"

#include <fun4all/Fun4AllServer.h>

#include <Geant4/G4ApplicationState.hh>  // for G4State_Idle
#include <Geant4/G4String.hh>
#include <Geant4/G4UIcmdWithAnInteger.hh>
#include <Geant4/G4UIdirectory.hh>  // for G4UIdirectory

#include <iosfwd>  // for std

class G4UIcommand;

using namespace std;

Fun4AllMessenger::Fun4AllMessenger(Fun4AllServer* ffa)
  : se(ffa)
{
  m_Fun4AllDir = new G4UIdirectory("/Fun4All/");
  m_Fun4AllDir->SetGuidance("UI commands to run Fun4All commands");
  m_RunCmd = new G4UIcmdWithAnInteger("/Fun4All/run", this);
  m_RunCmd->SetGuidance("Run Event(s)");
  m_RunCmd->SetParameterName("nEvents", 1);
  m_RunCmd->AvailableForStates(G4State_Idle);
}

Fun4AllMessenger::~Fun4AllMessenger()
{
  delete m_RunCmd;
  delete m_Fun4AllDir;
}

void Fun4AllMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == m_RunCmd)
  {
    se->run(m_RunCmd->GetNewIntValue(newValue));
  }
}
