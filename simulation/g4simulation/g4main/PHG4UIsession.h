// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4UISESSION_H
#define G4MAIN_PHG4UISESSION_H

#include <Geant4/G4UIsession.hh>

class PHG4UIsession : public G4UIsession {

public:
  PHG4UIsession();
  virtual ~PHG4UIsession() {}

  void Verbosity(int verb) {verbosity = verb;}
  
  G4UIsession * SessionStart();
  void PauseSessionStart(const G4String& Prompt);
  G4int ReceiveG4cout(const G4String& coutString);
  G4int ReceiveG4cerr(const G4String& cerrString);

private:
  int verbosity;

};

#endif //__PHG4UISESSION_H__
