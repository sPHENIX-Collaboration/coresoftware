#ifndef G4JET_CLUSTERJETINPUT_H
#define G4JET_CLUSTERJETINPUT_H

// first my own include
#include "JetInput.h"

// then other local incudes
#include "Jet.h"

// finally system includes
#include <vector>

// forward declarations
class PHCompositeNode;

class ClusterJetInput : public JetInput
{
 public:
  ClusterJetInput(Jet::SRC input);
  virtual ~ClusterJetInput() {}

  void identify(std::ostream& os = std::cout);

  Jet::SRC get_src() { return _input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode);

 private:
  int _verbosity;
  Jet::SRC _input;
};

#endif
