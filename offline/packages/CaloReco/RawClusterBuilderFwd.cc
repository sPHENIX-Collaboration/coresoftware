#include "RawClusterBuilderFwd.h"

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

class twrs_fwd
{
 public:
  explicit twrs_fwd(RawTower *);
  virtual ~twrs_fwd() = default;
  bool is_adjacent(const twrs_fwd &);
  void set_id(const int i)
  {
    id = i;
  }
  int get_id() const
  {
    return id;
  }
  int get_j_bin() const
  {
    return bin_j;
  }
  int get_k_bin() const
  {
    return bin_k;
  }

 protected:
  int bin_j;
  int bin_k;
  RawTowerDefs::keytype id;
};

twrs_fwd::twrs_fwd(RawTower *rt)
  : id(-1)
{
  bin_j = rt->get_bineta();
  bin_k = rt->get_binphi();
}

bool twrs_fwd::is_adjacent(twrs_fwd const &tower)
{
  if (bin_j - 1 <= tower.get_j_bin() && tower.get_j_bin() <= bin_j + 1)
  {
    if (bin_k - 1 <= tower.get_k_bin() && tower.get_k_bin() <= bin_k + 1)
    {
      return true;
    }
  }

  return false;
}

bool operator<(const twrs_fwd &a, const twrs_fwd &b)
{
  if (a.get_j_bin() != b.get_j_bin())
  {
    return a.get_j_bin() < b.get_j_bin();
  }
  return a.get_k_bin() < b.get_k_bin();
}

RawClusterBuilderFwd::RawClusterBuilderFwd(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , _min_tower_e(0.0)
  , chkenergyconservation(0)
  , detector("NONE")
{
}

int RawClusterBuilderFwd::InitRun(PHCompositeNode *topNode)
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

int RawClusterBuilderFwd::process_event(PHCompositeNode *topNode)
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
  std::vector<twrs_fwd> towerVector;
  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;
  for (; itr != begin_end.second; ++itr)
  {
    RawTower *tower = itr->second;
    RawTowerDefs::keytype towerid = itr->first;
    if (tower->get_energy() > _min_tower_e)
    {
      twrs_fwd twr(tower);
      twr.set_id(towerid);
      towerVector.push_back(twr);
    }
  }

  // cluster the towers
  std::multimap<int, twrs_fwd> clusteredTowers;
  PHMakeGroups(towerVector, clusteredTowers);

  RawCluster *cluster = nullptr;
  int last_id = -1;
  std::multimap<int, twrs_fwd>::iterator ctitr = clusteredTowers.begin();
  std::multimap<int, twrs_fwd>::iterator lastct = clusteredTowers.end();
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

    const twrs_fwd &tmptower = ctitr->second;
    RawTower *rawtower = towers->getTower(tmptower.get_id());

    const double e = rawtower->get_energy();
    cluster->addTower(rawtower->get_id(), e);
  }

  for (const auto &cluster_pair : _clusters->getClustersMap())
  {
    const RawClusterDefs::keytype clusterid = cluster_pair.first;
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
      std::cout << "RawClusterBuilder constucted ";
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

bool RawClusterBuilderFwd::CorrectPhi(RawCluster *cluster, RawTowerContainer *towers, RawTowerGeomContainer *towergeom)
{
  double sum = cluster->get_energy();
  double phimin = 999.;
  double phimax = -999.;
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  RawCluster::TowerConstIterator iter;
  for (iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    RawTower *tmpt = towers->getTower(iter->first);
    RawTowerGeom *tgeo =
        towergeom->get_tower_geometry(tmpt->get_id());
    double phi = tgeo->get_phi();
    if (phi > M_PI) phi = phi - 2. * M_PI;  // correct the cluster phi for slat geometry which is 0-2pi (L. Xue)
    if (phi < phimin)
    {
      phimin = phi;
    }
    if (phi > phimax)
    {
      phimax = phi;
    }
  }

  if ((phimax - phimin) < 3.) return false;  // cluster is not at phi discontinuity

  float mean = 0.;
  for (iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    RawTower *tmpt = towers->getTower(iter->first);
    double e = tmpt->get_energy();
    RawTowerGeom *tgeo =
        towergeom->get_tower_geometry(tmpt->get_id());
    double phi = tgeo->get_phi();
    if (phi > M_PI) phi = phi - 2. * M_PI;  // correct the cluster phi for slat geometry which is 0-2pi (L. Xue)
    if (phi < 0.)
    {
      phi = phi + 2. * M_PI;  // shift phi range for correct mean calculation
    }
    mean += e * phi;
  }
  mean = mean / sum;
  if (mean > M_PI)
  {
    mean = mean - 2. * M_PI;  // shift back
  }

  cluster->set_phi(mean);

  return true;  // mean phi was corrected
}

int RawClusterBuilderFwd::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderFwd::CreateNodes(PHCompositeNode *topNode)
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
