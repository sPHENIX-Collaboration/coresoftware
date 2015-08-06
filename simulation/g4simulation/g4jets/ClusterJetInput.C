
#include "ClusterJetInput.h"

#include "JetInput.h"
#include "Jet.h"
#include "JetV1.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

// standard includes
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

ClusterJetInput::ClusterJetInput(Jet::SRC input)
  : _verbosity(0),
    _input(input) {
}

std::vector<Jet*> ClusterJetInput::get_input(PHCompositeNode *topNode) {
  
  if (_verbosity > 0) cout << "ClusterJetInput::process_event -- entered" << endl;

  RawClusterContainer *clusters = NULL;
  if (_input == Jet::CEMC_CLUSTER) {
    clusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_CEMC");
    if (!clusters) {
      return std::vector<Jet*>();
    }
  } else if (_input == Jet::HCALIN_CLUSTER) {
    clusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALIN");
    if (!clusters) {
      return std::vector<Jet*>();
    }
  } else if (_input == Jet::HCALOUT_CLUSTER) {
    clusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALOUT");
    if (!clusters) {
      return std::vector<Jet*>();
    }
  } else {
    return std::vector<Jet*>();
  }
  
  std::vector<Jet*> pseudojets;
  RawClusterContainer::ConstRange begin_end = clusters->getClusters();
  RawClusterContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter !=  begin_end.second; ++rtiter) {
    RawCluster *cluster = rtiter->second;

    double eta = cluster->get_eta();
    double phi = cluster->get_phi();

    double pt = cluster->get_energy() / cosh(eta);
    double px = pt * cos(phi);
    double py = pt * sin(phi);
    double pz = pt * sinh(eta);

    Jet *jet = new JetV1();
    jet->set_px(px);
    jet->set_py(py);
    jet->set_pz(pz);
    jet->set_e(cluster->get_energy());
    jet->insert_comp(_input,0x0);
    pseudojets.push_back(jet);
  }

  if (_verbosity > 0) cout << "ClusterJetInput::process_event -- exited" << endl;

  return pseudojets;
}
