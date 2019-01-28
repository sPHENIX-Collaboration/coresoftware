#ifndef G4JET_TOWERJETINPUT_H
#define G4JET_TOWERJETINPUT_H

#include "JetInput.h"

#include "Jet.h"

#include <vector>

// forward declarations
class PHCompositeNode;

class TowerJetInput : public JetInput
{
 public:
  TowerJetInput(Jet::SRC input);
  virtual ~TowerJetInput() {}

  void identify(std::ostream& os = std::cout);

  Jet::SRC get_src() { return _input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode);

 private:
  int _verbosity;
  Jet::SRC _input;
};

#endif
