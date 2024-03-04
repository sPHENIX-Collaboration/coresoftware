#ifndef G4JET_JETPROBEINPUT__H
#define G4JET_JETPROBEINPUT__H

#include "JetInput.h"
#include "Jet.h"

#include <iostream>  // for cout, ostream
#include <vector>

#include <iostream>

class PHCompositeNode;

class JetProbeInput : public JetInput
{
 public:
  JetProbeInput(PHCompositeNode* node=nullptr);
  ~JetProbeInput() override {}

  //! by default, JetProbeInput process all truth primary particle.
  //! However, it can be configured to read only one or more embedded stream via add_embedding_flag()
  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return Jet::SRC::JET_PROBE; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

  // filled upon initialization if instantiated with topNode
  float phi { -100. };
  float eta { -100. };
  float pt  { -100. };

};

#endif
