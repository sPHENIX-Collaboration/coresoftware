#include "JetProbeInput.h"

#include "Jet.h"
#include "Jetv2.h"
#include "JetContainer.h"

#include <g4main/PHG4Particle.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// standard includes
#include <algorithm>  // std::find
#include <cmath>      // for asinh, sqrt
#include <cstdlib>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>

JetProbeInput::JetProbeInput(PHCompositeNode* topNode) {
  if (!topNode) { return; }

  JetContainer* jets = findNode::getClass<JetContainer>(topNode, "JetProbeContainer");
  if (!jets)
  {
    std::cout << PHWHERE << "JetProbeContainer node missing, doing nothing." << std::endl;
    return;
  }
  Jet* probe = jets->get_UncheckedAt(0);
  phi = probe->get_phi();
  eta = probe->get_eta();
  pt  = probe->get_pt();
}

void JetProbeInput::identify(std::ostream &os)
{
  os << "   JetProbeInput" << std::endl;
}

std::vector<Jet *> JetProbeInput::get_input(PHCompositeNode *topNode)
{
  if (Verbosity() > 0) { std::cout << "JetProbeInput::process_event -- entered" << std::endl; }

  JetContainer* jets = findNode::getClass<JetContainer>(topNode, "JetProbeContainer");
  if (!jets)
  {
    std::cout << PHWHERE << "JetProbeContainer node missing, doing nothing." << std::endl;
    return {};
  }

  // Pull the reconstructed track information off the node tree...
  std::vector<Jet *> pseudojets;
  Jet* jet = jets->get_UncheckedAt(0);
  Jet* probe = new Jetv2();
  probe->set_px(jet->get_px());
  probe->set_py(jet->get_py());
  probe->set_pz(jet->get_pz());
  probe->set_e(jet->get_e());
  pseudojets.push_back(probe);

  return pseudojets;
}
