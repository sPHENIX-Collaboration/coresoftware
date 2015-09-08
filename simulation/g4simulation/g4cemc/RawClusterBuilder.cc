#include "RawClusterBuilder.h"
#include "RawClusterContainer.h"
#include "RawClusterv1.h"
#include "PHMakeGroups.h"

#include "RawTower.h"
#include "RawTowerGeom.h"
#include "RawTowerContainer.h"

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>

// TODO: The output clusters assume the event is coming from (0,0,0) in the calculation of eta, phi...
// that should probably be modified at some point when things are more stable with the tracking vertex
// calculation.

using namespace std;

// this is just a helper class which enables us to handle rollovers
// when checking for adjacent towers, it requires one bit of
// information (the total number of phibins) which
// is not in the tower class
class twrs
{
public:
  twrs(RawTower *);
  virtual ~twrs() {}
  bool is_adjacent(twrs &);
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

twrs::twrs(RawTower *rt):
  maxphibin(-10),
  id(-1)
{
  bineta = rt->get_bineta();
  binphi = rt->get_binphi();
}

bool
twrs::is_adjacent(twrs &tower)
{
  if(bineta-1<=tower.get_bineta() && tower.get_bineta()<=bineta+1)
    {
      if(binphi-1<=tower.get_binphi() && tower.get_binphi()<=binphi+1)
	{
	  return true;
	}
      // cluster through the phi-wraparound
      else if(((tower.get_binphi() == maxphibin-1) && (binphi == 0)) ||
	      ((tower.get_binphi() == 0) && (binphi == maxphibin-1)))
	{
	  return true;
	}
    }
  
  return false;
}

bool operator<(const twrs& a, const twrs& b)
{
  if (a.get_bineta() != b.get_bineta())
    {
      return a.get_bineta() < b.get_bineta();
    }
  return a.get_binphi() < b.get_binphi();
}

RawClusterBuilder::RawClusterBuilder(const std::string& name):
  SubsysReco( name ),
  _clusters(NULL),
  _min_tower_e(0.0),
  chkenergyconservation(0),
  detector("NONE")
{}

int RawClusterBuilder::InitRun(PHCompositeNode *topNode)
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

int RawClusterBuilder::process_event(PHCompositeNode *topNode)
{

  string towernodename = "TOWER_" + detector;
  // Grab the towers
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, towernodename.c_str());
  if (!towers)
    {
      std::cout << PHWHERE << ": Could not find node " << towernodename.c_str() << std::endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeom *towergeom = findNode::getClass<RawTowerGeom>(topNode, towergeomnodename.c_str());
 if (! towergeom)
   {
     cout << PHWHERE << ": Could not find node " << towergeomnodename.c_str() << endl;
     return Fun4AllReturnCodes::ABORTEVENT;
   }
  // make the list of towers above threshold
  std::vector<twrs> towerVector;
  RawTowerContainer::ConstRange begin_end  = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;
  for (; itr != begin_end.second; ++itr)
    {
      RawTower* tower = itr->second;
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

  // extract the clusters
  std::vector<float> energy;
  std::vector<float> eta;
  std::vector<float> phi;

  std::multimap<int, twrs>::iterator ctitr = clusteredTowers.begin();
  std::multimap<int, twrs>::iterator lastct = clusteredTowers.end();
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

      twrs tmptower = ctitr->second;
      int iphi = tmptower.get_binphi();
      int ieta = tmptower.get_bineta();
      RawTower *rawtower = towers->getTower(ieta, iphi);
      if (tmptower.get_id() != (int) towers->genkey(ieta, iphi))
	{
	  cout << "id mismatch. internal: " << tmptower.get_id()
	       << ", towercontainer: " << towers->genkey(ieta, iphi)
	       << endl;
	  exit(1);
	}
      float e = rawtower->get_energy();
      energy[clusterid] += e;
      eta[clusterid] += e * towergeom->get_etacenter(rawtower->get_bineta());
      phi[clusterid] += e * towergeom->get_phicenter(rawtower->get_binphi());

      cluster->addTower(rawtower->get_id(), rawtower->get_energy());

      if (verbosity)
        {
          std::cout << "RawClusterBuilder id: " << (ctitr->first) << " Tower: "
                    << " (ieta,iphi) = (" << rawtower->get_bineta() << "," << rawtower->get_binphi() << ") "
                    << " (eta,phi,e) = (" << towergeom->get_etacenter(rawtower->get_bineta()) << ","
                    << towergeom->get_phicenter(rawtower->get_binphi()) << ","
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
          cout << "RawClusterBuilder: Cluster # " << icluster << " of " << nclusters << " "
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


bool RawClusterBuilder::CorrectPhi(RawCluster* cluster, RawTowerContainer* towers, RawTowerGeom *towergeom)
{
  double sum = cluster->get_energy();
  double phimin = 999.;
  double phimax = -999.;
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  RawCluster::TowerConstIterator iter;
  for (iter = begin_end.first; iter != begin_end.second; ++iter)
    { 
      RawTower* tmpt = towers->getTower(iter->first);
      double phi = towergeom->get_phicenter(tmpt->get_binphi());
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
      double phi = towergeom->get_phicenter(tmpt->get_binphi());
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


int RawClusterBuilder::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilder::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
    }

  _clusters = new RawClusterContainer();
  ClusterNodeName = "CLUSTER_" + detector;
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName.c_str(), "PHObject");
  dstNode->addNode(clusterNode);
}
