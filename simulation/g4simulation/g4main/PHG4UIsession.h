// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4UISESSION_H
#define G4MAIN_PHG4UISESSION_H

#include <Geant4/G4String.hh>     // for G4String
#include <Geant4/G4Types.hh>      // for G4int
#include <Geant4/G4UIsession.hh>

class PHG4UIsession : public G4UIsession {

public:
  PHG4UIsession();
  ~PHG4UIsession() override {}

  void Verbosity(int verb) {verbosity = verb;}
  
  G4UIsession * SessionStart() override;
  void PauseSessionStart(const G4String& Prompt) override;
  G4int ReceiveG4cout(const G4String& coutString) override;
  G4int ReceiveG4cerr(const G4String& cerrString) override;

private:
  int verbosity;

};

#endif //__PHG4UISESSION_H__
