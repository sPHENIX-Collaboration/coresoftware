#include "RawClusterBuilderFwd.h"
#include "RawClusterContainer.h"
#include "RawClusterv1.h"
#include "PHMakeGroups.h"

#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerGeom.h"
#include "RawTowerContainer.h"

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>

// TODO: The output clusters assume the event is coming from (0,0,0) in the calculation of eta, phi...
// that should probably be modified at some point when things are more stable with the tracking vertex
// calculation.

using namespace std;

class twrs_fwd
{
public:
  twrs_fwd(RawTower *);
  virtual ~twrs_fwd() {}
  bool is_adjacent(twrs_fwd &);
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

twrs_fwd::twrs_fwd(RawTower *rt):
  id(-1)
{
  bin_j = rt->get_bineta();
  bin_k = rt->get_binphi();
}

bool
twrs_fwd::is_adjacent(twrs_fwd &tower)
{
  if(bin_j-1<=tower.get_j_bin() && tower.get_j_bin()<=bin_j+1)
    {
      if(bin_k-1<=tower.get_k_bin() && tower.get_k_bin()<=bin_k+1)
	{
	  return true;
	}
    }
  
  return false;
}

bool operator<(const twrs_fwd& a, const twrs_fwd& b)
{
  if (a.get_j_bin() != b.get_j_bin())
    {
      return a.get_j_bin() < b.get_j_bin();
    }
  return a.get_k_bin() < b.get_k_bin();
}

RawClusterBuilderFwd::RawClusterBuilderFwd(const std::string& name):
  SubsysReco( name ),
  _clusters(NULL),
  _min_tower_e(0.0),
  chkenergyconservation(0),
  detector("NONE")
{}

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

  string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, towernodename.c_str());
  if (!towers)
    {
      std::cout << PHWHERE << ": Could not find node " << towernodename.c_str() << std::endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename.c_str());
  if (! towergeom)
   {
     cout << PHWHERE << ": Could not find node " << towergeomnodename.c_str() << endl;
     return Fun4AllReturnCodes::ABORTEVENT;
   }
  // make the list of towers above threshold
  std::vector<twrs_fwd> towerVector;
  RawTowerContainer::ConstRange begin_end  = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;
  for (; itr != begin_end.second; ++itr)
    {
      RawTower* tower = itr->second;
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

  // extract the clusters
  std::vector<float> energy;
  std::vector<float> eta;
  std::vector<float> phi;

  std::multimap<int, twrs_fwd>::iterator ctitr = clusteredTowers.begin();
  std::multimap<int, twrs_fwd>::iterator lastct = clusteredTowers.end();
  for (; ctitr != lastct; ++ctitr)
    {
      int clusterid = ctitr->first;
      RawCluster *cluster = _clusters->getCluster(clusterid);
      if (!cluster)
        {
          cluster = new RawClusterv1();
          _clusters->AddCluster(cluster);
          energy.push_back(0.0);
          eta.push_back(0.0);
          phi.push_back(0.0);
        }

      twrs_fwd tmptower = ctitr->second;
      int ij = tmptower.get_j_bin();
      int ik = tmptower.get_k_bin();
      RawTower *rawtower = towers->getTower(ij, ik);
      if (tmptower.get_id() != (int) RawTowerDefs::encode_towerid(towers->getCalorimeterID(), ij, ik ))
	{
	  cout <<__PRETTY_FUNCTION__<< " - Fatal Error - id mismatch. internal: " << tmptower.get_id()
	       << ", towercontainer: " << RawTowerDefs::encode_towerid( towers->getCalorimeterID(), ij, ik )
	       << endl;
	  exit(1);
	}
      float e = rawtower->get_energy();
      energy[clusterid] += e;
      RawTowerGeom *tgeo =  
	towergeom->get_tower_geometry(rawtower->get_id()); 
      eta[clusterid] += e * tgeo->get_eta();
      phi[clusterid] += e * tgeo->get_phi();

      cluster->addTower(rawtower->get_id(), rawtower->get_energy());

      if (verbosity)
        {
          std::cout << "RawClusterBuilderFwd id: " << (ctitr->first) << " Tower: "
                    << " (ieta,iphi) = (" << rawtower->get_bineta() << "," << rawtower->get_binphi() << ") "
                    << " (eta,phi,e) = (" << tgeo->get_eta() << ","
                    << tgeo->get_phi() << ","
                    << rawtower->get_energy() << ")"
                    << std::endl;
        }
    }

  unsigned int nclusters = _clusters->size();
  for (unsigned int icluster = 0; icluster < nclusters; icluster++)
    {
      if (energy[icluster] > 0)
        {
          eta[icluster] /= energy[icluster];
          phi[icluster] /= energy[icluster];
        }
      else
        {
          eta[icluster] = 0.0;
          phi[icluster] = 0.0;
        }

      if(phi[icluster] > M_PI)  phi[icluster] = phi[icluster] - 2.*M_PI; // convert [0,2pi] to [-pi,pi] for slat geometry(L. Xue)
      RawCluster *cluster = _clusters->getCluster(icluster);
      cluster->set_energy(energy[icluster]);
      cluster->set_eta(eta[icluster]);
      cluster->set_phi(phi[icluster]);

      if (verbosity)
        {
          cout << "RawClusterBuilderFwd: Cluster # " << icluster << " of " << nclusters << " "
                    << " (eta,phi,e) = (" << cluster->get_eta() << ", "
                    << cluster->get_phi() << ","
                    << cluster->get_energy() << ")"
                    << endl;
        }
    }

  // Correct mean Phi calculation for clusters at Phi discontinuity
  // Assumes that Phi goes from -pi to +pi
  for (unsigned int icluster = 0; icluster < nclusters; icluster++)
    {
      RawCluster *cluster = _clusters->getCluster(icluster);
      float oldphi = cluster->get_phi();
      bool corr = CorrectPhi(cluster, towers,towergeom);
      if (corr && verbosity)
        {
          std::cout << PHWHERE << " Cluster Phi corrected: " << oldphi << " " << cluster->get_phi() << std::endl;
        }
    }

  if (chkenergyconservation)
    {
      double ecluster = _clusters->getTotalEdep();
      double etower = towers->getTotalEdep();
      if (ecluster > 0)
	{
	  if (fabs(etower - ecluster) / ecluster > 1e-9)
	    {
	      cout << "energy conservation violation: ETower: " << etower
		   << " ECluster: " << ecluster 
		   << " diff: " << etower - ecluster << endl;
	    }
	}
      else
	{
	  if (etower != 0)
	    {
	      cout << "energy conservation violation: ETower: " << etower
                 << " ECluster: " << ecluster << endl;
	    }
	}
    }
  return Fun4AllReturnCodes::EVENT_OK;
}


bool RawClusterBuilderFwd::CorrectPhi(RawCluster* cluster, RawTowerContainer* towers, RawTowerGeomContainer *towergeom)
{
  double sum = cluster->get_energy();
  double phimin = 999.;
  double phimax = -999.;
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  RawCluster::TowerConstIterator iter;
  for (iter = begin_end.first; iter != begin_end.second; ++iter)
    { 
      RawTower* tmpt = towers->getTower(iter->first);
      RawTowerGeom *tgeo =  
	towergeom->get_tower_geometry(tmpt->get_id()); 
      double phi = tgeo->get_phi();
      if(phi > M_PI) phi = phi - 2.*M_PI; // correct the cluster phi for slat geometry which is 0-2pi (L. Xue)
      if (phi < phimin)
        {
          phimin = phi;
        }
      if (phi > phimax)
        {
          phimax = phi;
        }
    }

  if ((phimax - phimin) < 3.) return false; // cluster is not at phi discontinuity

  float mean = 0.;
  for (iter =begin_end.first; iter != begin_end.second; ++iter)
    { 
      RawTower* tmpt = towers->getTower(iter->first);
      double e = tmpt->get_energy();
      RawTowerGeom *tgeo =  
	towergeom->get_tower_geometry(tmpt->get_id()); 
      double phi = tgeo->get_phi();
      if(phi > M_PI) phi = phi - 2.*M_PI; // correct the cluster phi for slat geometry which is 0-2pi (L. Xue)
      if (phi < 0.)
        {
          phi = phi + 2.*M_PI;  // shift phi range for correct mean calculation
        }
      mean += e * phi;
    }
  mean = mean / sum;
  if (mean > M_PI)
    {
      mean = mean - 2.*M_PI;  // shift back
    }

  cluster->set_phi(mean);

  return true; // mean phi was corrected
}


int RawClusterBuilderFwd::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderFwd::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
    }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode",detector ));
  if(!DetNode){
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
   }

  _clusters = new RawClusterContainer();
  ClusterNodeName = "CLUSTER_" + detector;
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName.c_str(), "PHObject");
  DetNode->addNode(clusterNode);
}
