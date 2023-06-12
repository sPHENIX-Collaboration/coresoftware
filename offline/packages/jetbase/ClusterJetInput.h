#ifndef G4JET_CLUSTERJETINPUT_H
#define G4JET_CLUSTERJETINPUT_H

// first my own include
#include "JetInput.h"

// then other local incudes
#include "Jet.h"

// finally system includes
#include <iostream>  // for cout, ostream
#include <vector>

// forward declarations
class PHCompositeNode;

class ClusterJetInput : public JetInput
{
 public:
  ClusterJetInput(Jet::SRC input);
  ~ClusterJetInput() override {}

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return m_Input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

 private:
  int m_Verbosity = 0;
  Jet::SRC m_Input = Jet::VOID;
};

#endif
