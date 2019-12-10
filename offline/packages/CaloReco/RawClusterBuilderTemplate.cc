#include "RawClusterBuilderTemplate.h"

#include "BEmcCluster.h"
#include "BEmcRecCEMC.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>
#include <cstdio>
#include <exception>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

RawClusterBuilderTemplate::~RawClusterBuilderTemplate()
{
  delete bemc;
}

RawClusterBuilderTemplate::RawClusterBuilderTemplate(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , _min_tower_e(0.020)
  , chkenergyconservation(0)
  , detector("NONE")
{
  fEnergyNorm = 1.;
  bemc = new BEmcRecCEMC();
  //
  // Some initial values for clustering
  //
  // Configuration: number of towers in Phi and Eta, and tower dimension (here still in angle units and eta units)
  bemc->SetGeometry(262, 92, 1.0, 1.0);
  // Define vertex ... not used now
  float vertex[3] = {0, 0, 0};
  bemc->SetVertex(vertex);
  // Define threshold ... not used now
  bemc->SetTowerThreshold(0);
}

int RawClusterBuilderTemplate::InitRun(PHCompositeNode *topNode)
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

  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename.c_str());
  if (!towergeom)
  {
    cout << PHWHERE << ": Could not find node " << towergeomnodename.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int ixmin = 999999;
  int ixmax = -999999;
  int iymin = 999999;
  int iymax = -999999;
  RawTowerGeomContainer::ConstRange begin_end_geom = towergeom->get_tower_geometries();
  RawTowerGeomContainer::ConstIterator itr_geom = begin_end_geom.first;
  for (; itr_geom != begin_end_geom.second; ++itr_geom)
  {
    RawTowerGeom *towerg = itr_geom->second;
    RawTowerDefs::keytype towerid = towerg->get_id();
    int ix = RawTowerDefs::decode_index1(towerid);
    int iy = RawTowerDefs::decode_index2(towerid);
    if (ixmin > ix) ixmin = ix;
    if (ixmax < ix) ixmax = ix;
    if (iymin > iy) iymin = iy;
    if (iymax < iy) iymax = iy;
  }
  //  printf("************* Init CEMC: N of geom towers: %d; ix=%d-%d iy=%d-%d\n",ngeom,ixmin,ixmax,iymin,iymax);

  return Fun4AllReturnCodes::EVENT_OK;
}

bool RawClusterBuilderTemplate::Cell2Abs(RawTowerGeomContainer *towergeom, float phiC, float etaC, float &phi, float &eta)
{
  int NPHI = towergeom->get_phibins();
  int NETA = towergeom->get_etabins();

  int i1, i2;
  float dd;

  int iphi = phiC + 0.5;  // tower #
  if (iphi < 0 || iphi >= NPHI)
  {
    printf("RawClusterBuilderTemplate::Cell2Abs: wrong input phi: %d\n", iphi);
    return false;
  }

  i1 = iphi - 1;
  i2 = iphi + 1;
  if (i1 < 0)
    dd = fabs(towergeom->get_phicenter(i2) - towergeom->get_phicenter(iphi));
  else if (i2 >= NPHI)
    dd = fabs(towergeom->get_phicenter(iphi) - towergeom->get_phicenter(i1));
  else
    dd = fabs(towergeom->get_phicenter(i2) - towergeom->get_phicenter(i1)) / 2.;

  phi = towergeom->get_phicenter(iphi) + (phiC - iphi) * dd;

  int ieta = etaC + 0.5;  // tower #
  if (ieta < 0 || ieta >= NETA)
  {
    printf("RawClusterBuilderTemplate::Cell2Abs: wrong input eta: %d\n", ieta);
    return false;
  }

  i1 = ieta - 1;
  i2 = ieta + 1;
  if (i1 < 0)
    dd = fabs(towergeom->get_etacenter(i2) - towergeom->get_etacenter(ieta));
  else if (i2 >= NETA)
    dd = fabs(towergeom->get_etacenter(ieta) - towergeom->get_etacenter(i1));
  else
    dd = fabs(towergeom->get_etacenter(i2) - towergeom->get_etacenter(i1)) / 2.;

  eta = towergeom->get_etacenter(ieta) + (etaC - ieta) * dd;

  return true;
}

int RawClusterBuilderTemplate::process_event(PHCompositeNode *topNode)
{
  string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernodename.c_str());
  if (!towers)
  {
    std::cout << PHWHERE << ": Could not find node " << towernodename.c_str() << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename.c_str());
  if (!towergeom)
  {
    cout << PHWHERE << ": Could not find node " << towergeomnodename.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int NPHI = towergeom->get_phibins();
  int NETA = towergeom->get_etabins();
  /* 
 printf("Nphi= %d Neta= %d\n",NPHI,NETA);
 for( int i=0; i<NPHI; i++ ) printf("iphi= %3d %f\n",i,towergeom->get_phicenter(i));
 for( int i=0; i<NETA; i++ ) {
   //   std::pair<double, double> beta0 = towergeom->get_etabounds(i);
   printf("ieta= %3d %f\n",i,towergeom->get_etacenter(i));
   //   printf("ieta= %3d %f %f %f\n",i,towergeom->get_etacenter(i),beta0.first,beta0.second);
 }
 */

  bemc->SetGeometry(NPHI, NETA, 2. * M_PI / NPHI, 0.024);  // !!!!! The last parameter not used for now

  // _clusters->Reset(); // !!! Not sure if it is necessarry to do it - ask Chris

  // make the list of towers above threshold
  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;

  // Define vector of towers in EmcModule format to input into BEmc
  EmcModule vhit;
  std::vector<EmcModule> HitList;
  HitList.erase(HitList.begin(), HitList.end());
  int ich, ix, iy;

  for (; itr != begin_end.second; ++itr)
  {
    RawTower *tower = itr->second;
    //      printf("  Tower e=%f (%f)\n",tower->get_energy(), _min_tower_e);
    if (tower->get_energy() > _min_tower_e)
    {
      ix = tower->get_binphi();
      iy = tower->get_bineta();
      ich = iy * NPHI + ix;
      vhit.ich = ich;
      vhit.amp = tower->get_energy() * fEnergyNorm;  // !!! Global Calibration
      vhit.tof = 0.;
      HitList.push_back(vhit);
    }
  }

  bemc->SetModules(&HitList);

  // Find clusters (as a set of towers with common edge)
  int ncl = bemc->FindClusters();
  if (ncl < 0)
  {
    printf("!!! Error in BEmcRec::FindClusters(): Too many clusters, fgMaxLen parameter needs to be increased\n");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get pointer to clusters
  std::vector<EmcCluster> *ClusterList = bemc->GetClusters();
  std::vector<EmcCluster>::iterator pc;

  std::vector<EmcPeakarea>::iterator pp;
  float ecl, ecore, xcg, ycg, xx, xy, yy;
  float xcorr, ycorr;
  EmcModule hmax;
  RawCluster *cluster;

  std::vector<EmcPeakarea> PList;
  std::vector<EmcModule> Peaks;
  std::vector<EmcPeakarea> *pPList = &PList;
  std::vector<EmcModule> *pPeaks = &Peaks;

  int iphi, ieta;
  float phi, eta;
  float prob, chi2;
  int ndf;

  vector<EmcModule>::iterator ph;
  vector<EmcModule> hlist;

  // ncl = 0;
  for (pc = ClusterList->begin(); pc != ClusterList->end(); ++pc)
  {
    //    ecl = pc->GetTotalEnergy();
    //    pc->GetMoments( &xcg, &ycg, &xx, &xy, &yy );

    int npk = pc->GetPeaks(pPList, pPeaks);
    if (npk < 0) return Fun4AllReturnCodes::ABORTEVENT;

    //    printf("  iCl=%d (%d): E=%f  x=%f  y=%f\n",ncl,npk,ecl,xcg,ycg);

    for (pp = pPList->begin(); pp != pPList->end(); ++pp)
    {
      // Cluster energy
      ecl = pp->GetTotalEnergy();
      ecore = pp->GetECoreCorrected();
      // 3x3 energy around center of gravity
      //e9 = pp->GetE9();
      // Ecore (basically near 2x2 energy around center of gravity)
      //ecore = pp->GetECore();
      // Center of Gravity etc.
      pp->GetMoments(&xcg, &ycg, &xx, &xy, &yy);
      // Tower with max energy
      hmax = pp->GetMaxTower();

      //      phi = (xcg-float(NPHI)/2.+0.5)/float(NPHI)*2.*M_PI;
      //      eta = (ycg-float(NETA)/2.+0.5)/float(NETA)*2.2; // -1.1<eta<1.1;

      //      Cell2Abs(towergeom,xcg,ycg,phi,eta);

      pp->GetCorrPos(&xcorr, &ycorr);
      Cell2Abs(towergeom, xcorr, ycorr, phi, eta);
      const double ref_radius = towergeom->get_radius();

      if (phi > M_PI) phi -= 2. * M_PI;  // convert to [-pi,pi]]

      prob = pp->GetProb(chi2, ndf);
      //      printf("Prob/Chi2/NDF= %f %f %d Ecl=%f\n",prob,chi2,ndf,ecl);

      cluster = new RawClusterv1();
      cluster->set_energy(ecl);
      cluster->set_ecore(ecore);
      cluster->set_phi(phi);
      cluster->set_r(ref_radius);
      cluster->set_z(ref_radius * sinh(eta));
      cluster->set_prob(prob);
      if (ndf > 0)
        cluster->set_chi2(chi2 / ndf);
      else
        cluster->set_chi2(0);

      hlist = pp->GetHitList();
      ph = hlist.begin();
      while (ph != hlist.end())
      {
        ich = (*ph).ich;
        ieta = ich / NPHI;
        iphi = ich % NPHI;
        // that code needs a closer look - here are the towers
        // with their energy added to the cluster object where
        // the id is the tower id
        // !!!!! Make sure twrkey is correctly extracted
        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(towers->getCalorimeterID(), ieta, iphi);
        //	printf("%d %d: %d e=%f\n",iphi,ieta,twrkey,(*ph).amp);
        cluster->addTower(twrkey, (*ph).amp / fEnergyNorm);
        ++ph;
      }

      _clusters->AddCluster(cluster);
      // ncl++;

      //      printf("    ipk=%d: E=%f  E9=%f  x=%f  y=%f  MaxTower: (%d,%d) e=%f\n",ipk,ecl,e9,xcg,ycg,hmax.ich%NPHI,hmax.ich/NPHI,hmax.amp);
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

bool RawClusterBuilderTemplate::CorrectPhi(RawCluster *cluster, RawTowerContainer *towers, RawTowerGeomContainer *towergeom)
{
  double sum = cluster->get_energy();
  double phimin = 999.;
  double phimax = -999.;
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  RawCluster::TowerConstIterator iter;
  for (iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    RawTower *tmpt = towers->getTower(iter->first);
    double phi = towergeom->get_phicenter(tmpt->get_binphi());
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
    double phi = towergeom->get_phicenter(tmpt->get_binphi());
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

int RawClusterBuilderTemplate::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTemplate::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  //Get the _det_name subnode
  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", detector));

  //Check that it is there
  if (!cemcNode)
  {
    cemcNode = new PHCompositeNode(detector);
    dstNode->addNode(cemcNode);
  }

  _clusters = new RawClusterContainer();
  ClusterNodeName = "CLUSTER_" + detector;
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName.c_str(), "PHObject");
  cemcNode->addNode(clusterNode);
}
