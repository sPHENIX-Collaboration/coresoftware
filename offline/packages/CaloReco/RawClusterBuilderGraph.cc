#include "RawClusterBuilderGraph.h"

#include "PHMakeGroups.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

// this is just a helper class which enables us to handle rollovers
// when checking for adjacent towers, it requires one bit of
// information (the total number of phibins) which
// is not in the tower class
class twrs
{
 public:
  explicit twrs(RawTower *);
  virtual ~twrs() = default;
  bool is_adjacent(const twrs &);
  void set_id(const int i)
  {
    id = i;
  }
  int get_id() const
  {
    return id;
  }
  void set_maxphibin(const int i)
  {
    maxphibin = i;
  }
  int get_maxphibin() const
  {
    return maxphibin;
  }
  int get_bineta() const
  {
    return bineta;
  }
  int get_binphi() const
  {
    return binphi;
  }

 protected:
  int bineta;
  int binphi;
  int maxphibin;
  RawTowerDefs::keytype id;
};

twrs::twrs(RawTower *rt)
  : maxphibin(-10)
  , id(-1)
{
  bineta = rt->get_bineta();
  binphi = rt->get_binphi();
}

bool twrs::is_adjacent(const twrs &tower)
{
  if (bineta - 1 <= tower.get_bineta() && tower.get_bineta() <= bineta + 1)
  {
    if (binphi - 1 <= tower.get_binphi() && tower.get_binphi() <= binphi + 1)
    {
      return true;
    }
    // cluster through the phi-wraparound
    else if (((tower.get_binphi() == maxphibin - 1) && (binphi == 0)) ||
             ((tower.get_binphi() == 0) && (binphi == maxphibin - 1)))
    {
      return true;
    }
  }

  return false;
}

bool operator<(const twrs &a, const twrs &b)
{
  if (a.get_bineta() != b.get_bineta())
  {
    return a.get_bineta() < b.get_bineta();
  }
  return a.get_binphi() < b.get_binphi();
}

RawClusterBuilderGraph::RawClusterBuilderGraph(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , _min_tower_e(0.0)
  , chkenergyconservation(0)
  , detector("NONE")
{
}

int RawClusterBuilderGraph::InitRun(PHCompositeNode *topNode)
{
  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderGraph::process_event(PHCompositeNode *topNode)
{
  std::string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernodename);
  if (!towers)
  {
    std::cout << PHWHERE << ": Could not find node " << towernodename << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  std::string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnodename << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  // make the list of towers above threshold
  std::vector<twrs> towerVector;
  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;
  for (; itr != begin_end.second; ++itr)
  {
    RawTower *tower = itr->second;
    RawTowerDefs::keytype towerid = itr->first;
    if (tower->get_energy() > _min_tower_e)
    {
      twrs twr(tower);
      twr.set_maxphibin(towergeom->get_phibins());
      twr.set_id(towerid);
      towerVector.push_back(twr);
    }
  }

  // cluster the towers
  std::multimap<int, twrs> clusteredTowers;
  PHMakeGroups(towerVector, clusteredTowers);

  RawCluster *cluster = nullptr;
  int last_id = -1;
  std::multimap<int, twrs>::iterator ctitr = clusteredTowers.begin();
  std::multimap<int, twrs>::iterator lastct = clusteredTowers.end();
  for (; ctitr != lastct; ++ctitr)
  {
    int clusterid = ctitr->first;

    if (last_id != clusterid)
    {
      // new cluster
      cluster = new RawClusterv1();
      _clusters->AddCluster(cluster);

      last_id = clusterid;
    }
    assert(cluster);

    const twrs &tmptower = ctitr->second;
    RawTower *rawtower = towers->getTower(tmptower.get_id());

    const double e = rawtower->get_energy();
    cluster->addTower(rawtower->get_id(), e);
  }

  for (const auto &cluster_pair : _clusters->getClustersMap())
  {
    RawClusterDefs::keytype clusterid = cluster_pair.first;
    RawCluster *clusterA = cluster_pair.second;

    assert(clusterA);
    assert(clusterA->get_id() == clusterid);

    double sum_x(0);
    double sum_y(0);
    double sum_z(0);
    double sum_e(0);

    for (const auto tower_pair : clusterA->get_towermap())
    {
      const RawTower *rawtower = towers->getTower(tower_pair.first);
      const RawTowerGeom *rawtowergeom = towergeom->get_tower_geometry(tower_pair.first);

      assert(rawtower);
      assert(rawtowergeom);
      const double e = rawtower->get_energy();

      sum_e += e;

      if (e > 0)
      {
        sum_x += e * rawtowergeom->get_center_x();
        sum_y += e * rawtowergeom->get_center_y();
        sum_z += e * rawtowergeom->get_center_z();
      }
    }  //     for (const auto tower_pair : clusterA->get_towermap())

    clusterA->set_energy(sum_e);

    if (sum_e > 0)
    {
      sum_x /= sum_e;
      sum_y /= sum_e;
      sum_z /= sum_e;

      clusterA->set_r(sqrt(sum_y * sum_y + sum_x * sum_x));
      clusterA->set_phi(atan2(sum_y, sum_x));
      clusterA->set_z(sum_z);
    }

    if (Verbosity() > 1)
    {
      std::cout << "RawClusterBuilderGraph constucted ";
      clusterA->identify();
    }
  }  //  for (const auto & cluster_pair : _clusters->getClustersMap())

  if (chkenergyconservation)
  {
    double ecluster = _clusters->getTotalEdep();
    double etower = towers->getTotalEdep();
    if (ecluster > 0)
    {
      if (fabs(etower - ecluster) / ecluster > 1e-9)
      {
        std::cout << "energy conservation violation: ETower: " << etower
                  << " ECluster: " << ecluster
                  << " diff: " << etower - ecluster << std::endl;
      }
    }
    else
    {
      if (etower != 0)
      {
        std::cout << "energy conservation violation: ETower: " << etower
                  << " ECluster: " << ecluster << std::endl;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderGraph::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderGraph::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(detector);
    dstNode->addNode(DetNode);
  }

  _clusters = new RawClusterContainer();
  ClusterNodeName = "CLUSTER_" + detector;
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName, "PHObject");
  DetNode->addNode(clusterNode);
}
