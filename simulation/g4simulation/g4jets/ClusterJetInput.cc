
#include "ClusterJetInput.h"

#include "Jet.h"
#include "Jetv1.h"

#include <phool/getClass.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector

// standard includes
#include <cassert>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>

using namespace std;

ClusterJetInput::ClusterJetInput(Jet::SRC input)
  : _verbosity(0)
  , _input(input)
{
}

void ClusterJetInput::identify(std::ostream &os)
{
  os << "   ClusterJetInput: ";
  if (_input == Jet::CEMC_CLUSTER)
    os << "CLUSTER_CEMC to Jet::CEMC_CLUSTER";
  if (_input == Jet::EEMC_CLUSTER)
    os << "CLUSTER_EEMC to Jet::EEMC_CLUSTER";
  else if (_input == Jet::HCALIN_CLUSTER)
    os << "CLUSTER_HCALIN to Jet::HCALIN_CLUSTER";
  else if (_input == Jet::HCALOUT_CLUSTER)
    os << "CLUSTER_HCALOUT to Jet::HCALOUT_CLUSTER";
  os << endl;
}

std::vector<Jet *> ClusterJetInput::get_input(PHCompositeNode *topNode)
{
  if (_verbosity > 0) cout << "ClusterJetInput::process_event -- entered" << endl;
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    cout << "ClusterJetInput::get_input - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    assert(vertexmap);  // force quit

    return std::vector<Jet *>();
  }

  if (vertexmap->empty())
  {
    cout << "ClusterJetInput::get_input - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_bbc or tracking reco flags in the main macro in order to reconstruct the global vertex." << endl;
    return std::vector<Jet *>();
  }

  RawClusterContainer *clusters = nullptr;
  if (_input == Jet::CEMC_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::EEMC_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_EEMC");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::HCALIN_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::HCALOUT_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::HCAL_TOPO_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_HCAL");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::ECAL_TOPO_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_EMCAL");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::FEMC_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_FEMC");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else if (_input == Jet::FHCAL_CLUSTER)
  {
    clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_FHCAL");
    if (!clusters)
    {
      return std::vector<Jet *>();
    }
  }
  else
  {
    return std::vector<Jet *>();
  }

  // first grab the event vertex or bail
  GlobalVertex *vtx = vertexmap->begin()->second;
  CLHEP::Hep3Vector vertex;
  if (vtx)
    vertex.set(vtx->get_x(), vtx->get_y(), vtx->get_z());
  else
    return std::vector<Jet *>();

  std::vector<Jet *> pseudojets;
  RawClusterContainer::ConstRange begin_end = clusters->getClusters();
  RawClusterContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    RawCluster *cluster = rtiter->second;

    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*cluster, vertex);

    Jet *jet = new Jetv1();
    jet->set_px(E_vec_cluster.x());
    jet->set_py(E_vec_cluster.y());
    jet->set_pz(E_vec_cluster.z());
    jet->set_e(cluster->get_energy());
    jet->insert_comp(_input, cluster->get_id());
    pseudojets.push_back(jet);
  }

  if (_verbosity > 0) cout << "ClusterJetInput::process_event -- exited" << endl;

  return pseudojets;
}
