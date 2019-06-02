#include "PHG4UIsession.h"
 
#include <iostream>  // for flush, basic_ostream, ostream, cerr, cout
#include <string>    // for operator<<

PHG4UIsession::PHG4UIsession()
  : verbosity(0) {
}

G4UIsession * PHG4UIsession::SessionStart() { return nullptr; }
 
void PHG4UIsession::PauseSessionStart(const G4String&) {;}

G4int PHG4UIsession::ReceiveG4cout(const G4String& coutString) {
  //  if (verbosity > 0) std::cout << coutString << std::flush;
  std::cout << coutString << std::flush;
  return 0;
}
 
G4int PHG4UIsession::ReceiveG4cerr(const G4String& cerrString) {
  std::cerr << cerrString << std::flush;
  return 0;
}                              
