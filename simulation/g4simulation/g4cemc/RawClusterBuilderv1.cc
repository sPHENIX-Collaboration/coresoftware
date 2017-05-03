#include "RawClusterBuilderv1.h"
#include "RawClusterContainer.h"
#include "RawClusterv1.h"
#include "PHMakeGroups.h"

#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerContainer.h"

#include "BEmcRec.h"
#include "BEmcCluster.h"

#define MaxNofPeaks 100

#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>

using namespace std;


RawClusterBuilderv1::~RawClusterBuilderv1()
{
  delete bemc;
}

RawClusterBuilderv1::RawClusterBuilderv1(const std::string& name):
  SubsysReco( name ),
  _clusters(NULL),
  _min_tower_e(0.0),
  chkenergyconservation(0),
  detector("NONE")
{
  fEnergyNorm = 1.;
  bemc = new BEmcRec();
  //
  // Some initial values for clustering
  //
  // Configuration: number of towers in Phi and Eta, and tower dimension (here still in angle units and eta units)
  bemc->SetGeometry( 262, 92, 1.0, 1.0 );
  // Define vertex ... not used now
  float vertex[3] = {0,0,0};
  bemc->SetVertex(vertex);
  // Define threshold ... not used now
  bemc->SetTowerThreshold(0);
}

int RawClusterBuilderv1::InitRun(PHCompositeNode *topNode)
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

int RawClusterBuilderv1::process_event(PHCompositeNode *topNode)
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

 int NPHI = towergeom->get_phibins();
 int NETA = towergeom->get_etabins();
 bemc->SetConf( NPHI, NETA );
 _clusters->Reset(); // !!! Not sure if it is necessarry to do it - ask Chris

  // make the list of towers above threshold
  RawTowerContainer::ConstRange begin_end  = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;

  // Define vector of towers in EmcModule format to input into BEmc
  EmcModule vhit;
  std::vector<EmcModule> HitList;
  HitList.erase(HitList.begin(), HitList.end());
  int ich, ix, iy;

  for (; itr != begin_end.second; ++itr)
    {
      RawTower* tower = itr->second;
      if (tower->get_energy() > _min_tower_e)
        {
	  ix = tower->get_binphi();
	  iy = tower->get_bineta();
	  ich = iy*NPHI + ix;
	  vhit.ich = ich;
	  vhit.amp = tower->get_energy()*fEnergyNorm; // !!! Global Calibration
	  vhit.tof = 0.;
	  HitList.push_back(vhit);
        }
    }

  bemc->SetModules(&HitList);

  // Get pointer to clusters
  std::vector<EmcCluster> *ClusterList = bemc->GetClusters();
  std::vector<EmcCluster>::iterator pc;

  float ecl, xcg, ycg, xx, xy, yy;
  EmcPeakarea pPList[MaxNofPeaks];
  EmcPeakarea *pp;
  EmcModule peaks[MaxNofPeaks];
  EmcModule hmax;
  RawCluster *cluster;

  int iphi, ieta;
  float phi, eta;
  float dphi, phistep, deta, etastep;

  vector<EmcModule>::iterator ph;
  vector<EmcModule> hlist;

  int ncl = 0;
  for( pc=ClusterList->begin(); pc!=ClusterList->end(); ++pc){

    //    ecl = pc->GetTotalEnergy();
    //    pc->GetMoments( &xcg, &ycg, &xx, &xy, &yy );
    //    printf("Cl: %f %f\n",xcg,ycg);
    int npk = pc->GetPeaks(pPList, peaks);
    pp = pPList;

    //    printf("  iCl=%d  NPeaks=%d: E=%f  x=%f  y=%f\n",iCl,npk,ecl,xcg,ycg);

    for(int ipk=0; ipk<npk; ipk++){

      // Cluster energy
      ecl = pp->GetTotalEnergy();
      // 3x3 energy around center of gravity
      //e9 = pp->GetE9();
      // Ecore (basically near 2x2 energy around center of gravity)
      //ecore = pp->GetECore();
      // Center of Gravity etc.
      pp->GetMoments( &xcg, &ycg, &xx, &xy, &yy );
      // Tower with max energy
      hmax = pp->GetMaxTower();

      //      phi = (xcg-float(NPHI)/2.+0.5)/float(NPHI)*2.*M_PI;
      //      eta = (ycg-float(NETA)/2.+0.5)/float(NETA)*2.2; // -1.1<eta<1.1;

      iphi = xcg+0.5;
      dphi = xcg - float(iphi); // this is from -0.5 to +0.5
      phi = towergeom->get_phicenter(iphi);
      std::pair<double, double> phibounds = towergeom->get_phibounds(iphi);
      phistep = phibounds.second - phibounds.first;
      phi += dphi*phistep;

      ieta = ycg+0.5;
      deta = ycg - float(ieta); // this is from -0.5 to +0.5
      eta = towergeom->get_etacenter(ieta);
      //      etastep = towergeom->get_etastep();
      std::pair<double, double> beta = towergeom->get_etabounds(ieta);
      etastep = fabs(beta.second-beta.first);
      eta += deta*etastep;

      cluster = new RawClusterv1();
      cluster->set_energy(ecl);
      cluster->set_phi(phi);
      cluster->set_eta(eta);

      hlist = pp->GetHitList();
      ph = hlist.begin();
      while( ph != hlist.end() ) {
	ich = (*ph).ich;
	ieta = ich/NPHI;
	iphi = ich%NPHI;
	// that code needs a closer look - here are the towers 
	// with their energy added to the cluster object where 
	// the id is the tower id
	RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid( RawTowerDefs::NONE , ieta , iphi );
	cluster->addTower(twrkey,(*ph).amp/fEnergyNorm);
	++ph;
      }

      _clusters->AddCluster(cluster);
      ncl++;

      //      printf("    ipk=%d: E=%f  E9=%f  x=%f  y=%f  MaxTower: (%d,%d) e=%f\n",ipk,ecl,e9,xcg,ycg,hmax.ich%NPHI,hmax.ich/NPHI,hmax.amp);
      pp++;

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


bool RawClusterBuilderv1::CorrectPhi(RawCluster* cluster, RawTowerContainer* towers, RawTowerGeomContainer *towergeom)
{
  double sum = cluster->get_energy();
  double phimin = 999.;
  double phimax = -999.;
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  RawCluster::TowerConstIterator iter;
  for (iter =begin_end.first; iter != begin_end.second; ++iter)
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


int RawClusterBuilderv1::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderv1::CreateNodes(PHCompositeNode *topNode)
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
