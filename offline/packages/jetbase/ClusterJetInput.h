#ifndef JETBASE_CLUSTERJETINPUT_H
#define JETBASE_CLUSTERJETINPUT_H

// first my own include
#include "JetInput.h"

// then other local incudes
#include "Jet.h"

#include <globalvertex/GlobalVertex.h>

// finally system includes
#include <iostream>  // for cout, ostream
#include <vector>

// forward declarations
class PHCompositeNode;
class GlobalVertex;

class ClusterJetInput : public JetInput
{
 public:
  ClusterJetInput(Jet::SRC input);
  ~ClusterJetInput() override {}

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return m_Input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

  void set_GlobalVertexType(GlobalVertex::VTXTYPE type) 
  {
    m_use_vertextype = true;
    m_vertex_type = type;
  }

 private:
  int m_Verbosity = 0;
  Jet::SRC m_Input = Jet::VOID;
  bool m_use_vertextype {false};
  GlobalVertex::VTXTYPE m_vertex_type = GlobalVertex::UNDEFINED;
};

#endif
