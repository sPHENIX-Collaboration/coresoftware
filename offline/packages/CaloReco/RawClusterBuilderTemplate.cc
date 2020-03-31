#include "RawClusterBuilderTemplate.h"

#include "BEmcCluster.h"
#include "BEmcRec.h"
#include "BEmcRecFEMC.h"
#include "BEmcRecCEMC.h"
#include "BEmcRecEEMC.h"
#include "BEmcProfile.h"

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
  if( bemc != nullptr ) delete bemc;
}

RawClusterBuilderTemplate::RawClusterBuilderTemplate(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , _min_tower_e(0.020)
  , chkenergyconservation(0)
  , detector("NONE")
  , bPrintGeom(false)
{
  bemc = nullptr;

  BINX0 = 0;
  NBINX = 0;
  BINY0 = 0;
  NBINY = 0;

  fEnergyNorm = 1.;
}

void RawClusterBuilderTemplate::Detector(const std::string& d)
{ 
  detector = d;

  // Create proper BEmcRec object

  if( detector == "CEMC" ) {
    bemc = new BEmcRecCEMC();
  }
  else if( detector == "FEMC" ) {
    bemc = new BEmcRecFEMC();
  }
  else if( detector == "EEMC" ) {
    bemc = new BEmcRecEEMC();
  }
  else {
    printf("Warning from RawClusterBuilderTemplate::Detector(): no detector specific class BEmcRec{Name} defined for detector %s. Default BEmcRec will be used\n",detector.c_str());
    bemc = new BEmcRec();
  }

  // Define vertex ... not used now
  float vertex[3] = {0, 0, 0};
  bemc->SetVertex(vertex);
  // Define threshold ... not used now
  bemc->SetTowerThreshold(0);
}

void RawClusterBuilderTemplate::LoadProfile(const char *fname) 
{ 
  if( bemc == nullptr ) {
    printf("Error in RawClusterBuilderTemplate::LoadProfile(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it\n");
    return;
  }

  //  _emcprof = new BEmcProfile(fname); 
  bemc->LoadProfile(fname);
}

void RawClusterBuilderTemplate::SetCylindricalGeometry()
{ 
  if( bemc == nullptr ) {
    printf("Error in RawClusterBuilderTemplate::SetCylindricalGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it\n");
    return;
  }

  bemc->SetCylindricalGeometry();
}

void RawClusterBuilderTemplate::SetPlanarGeometry()
{ 
  if( bemc == nullptr ) {
    printf("Error in RawClusterBuilderTemplate::SetPlanarGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it\n");
    return;
  }

  bemc->SetPlanarGeometry();
}

int RawClusterBuilderTemplate::InitRun(PHCompositeNode *topNode)
{
  if( bemc == nullptr ) {
    printf("Error in RawClusterBuilderTemplate::InitRun(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it\n");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

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

  int ngeom = 0;
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
    int ix = RawTowerDefs::decode_index2(towerid); // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(towerid); // index1 is eta in CYL
    if (ixmin > ix) ixmin = ix;
    if (ixmax < ix) ixmax = ix;
    if (iymin > iy) iymin = iy;
    if (iymax < iy) iymax = iy;
    ngeom++;
  }
  printf("Infor from RawClusterBuilderTemplate::InitRun(): Init geometry for %s: N of geom towers: %d; ix=%d-%d iy=%d-%d\n",detector.c_str(),ngeom,ixmin,ixmax,iymin,iymax);

  if (ixmax < ixmin || iymax < iymin) {
    printf("Error in RawClusterBuilderTemplate::InitRun(): wrong geometry data for detector %s\n",detector.c_str());
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  BINX0 = ixmin;
  NBINX = ixmax - ixmin + 1;
  BINY0 = iymin;
  NBINY = iymax - iymin + 1;

  bemc->SetGeometry(NBINX, NBINY, 1, 1);  // !!!!! The last parameter not used for now

  itr_geom = begin_end_geom.first;
  for (; itr_geom != begin_end_geom.second; ++itr_geom)
  {
    RawTowerGeom *towerg = itr_geom->second;
    RawTowerDefs::keytype towerid = towerg->get_id();
    //    int itype = towerg->get_tower_type();
    //    if( itype==2 ) { // PbSc
    int ix = RawTowerDefs::decode_index2(towerid); // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(towerid); // index1 is eta in CYL
    ix -= BINX0;
    iy -= BINY0;
    bemc->SetTowerGeometry(ix, iy, towerg->get_center_x(), towerg->get_center_y(), towerg->get_center_z());
    //    }
  }

  if(!bemc->CompleteTowerGeometry()) return Fun4AllReturnCodes::ABORTEVENT;

  if(bPrintGeom) {
    string fname = "geom_" + detector + ".txt";
    //    bemc->PrintTowerGeometry("geom.txt");
    bemc->PrintTowerGeometry(fname.c_str());
    //    PrintCylGeom(towergeom,"phieta.txt");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTemplate::PrintCylGeom(RawTowerGeomContainer *towergeom, const char* fname)
{
  FILE* pf = fopen(fname, "w");
  if (!pf)
  {
    printf("Error in BEmcRec::RawClusterBuilderTemplate::PrintCylGeom(): Failed to open file %s\n", fname);
    return;
  }

  fprintf(pf, "%d %d\n", NBINX, NBINY);
  for( int ip=0; ip<NBINX; ip++ ) {
    float phi = towergeom->get_phicenter(ip);
    fprintf(pf, "%d %f\n", ip, phi);
  }
  for( int ip=0; ip<NBINY; ip++ ) {
    float eta = towergeom->get_etacenter(ip);
    fprintf(pf, "%d %f\n", ip, eta);
  }

  fclose(pf);
}


bool RawClusterBuilderTemplate::Cell2Abs(RawTowerGeomContainer *towergeom, float phiC, float etaC, float &phi, float &eta)
{
  phi = eta = 0;
  return false;
}

int RawClusterBuilderTemplate::process_event(PHCompositeNode *topNode)
{
  if( bemc == nullptr ) {
    printf("Error in RawClusterBuilderTemplate::process_event(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it\n");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

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
      //	  printf("(%d,%d)  (%d,%d)\n",tower->get_column(),tower->get_row(),tower->get_binphi(),tower->get_bineta());
      //	  ix = tower->get_column();
      RawTowerDefs::keytype towerid = tower->get_id();
      int ix = RawTowerDefs::decode_index2(towerid); // index2 is phi in CYL
      int iy = RawTowerDefs::decode_index1(towerid); // index1 is eta in CYL
      ix -= BINX0; 
      iy -= BINY0;
      //      ix = tower->get_bineta() - BINX0;  // eta: index1
      //      iy = tower->get_binphi() - BINY0;  // phi: index2
      if (ix >= 0 && ix < NBINX && iy >= 0 && iy < NBINY)
      {
        ich = iy * NBINX + ix;
        vhit.ich = ich;
        vhit.amp = tower->get_energy() * fEnergyNorm;  // !!! Global Calibration
        vhit.tof = 0.;
        HitList.push_back(vhit);
      }
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
  //  float xcorr, ycorr;
  EmcModule hmax;
  RawCluster *cluster;

  std::vector<EmcPeakarea> PList;
  std::vector<EmcModule> Peaks;
  std::vector<EmcPeakarea> *pPList = &PList;
  std::vector<EmcModule> *pPeaks = &Peaks;

  float prob, chi2;
  int ndf;
  float xg, yg, zg;

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
      pp->GetGlobalPos(xg, yg, zg);

      // Tower with max energy
      hmax = pp->GetMaxTower();

      //      phi = (xcg-float(NPHI)/2.+0.5)/float(NPHI)*2.*M_PI;
      //      eta = (ycg-float(NETA)/2.+0.5)/float(NETA)*2.2; // -1.1<eta<1.1;

      //      Cell2Abs(towergeom,xcg,ycg,phi,eta);

      //      pp->GetCorrPos(&xcorr, &ycorr);
      //      Cell2Abs(towergeom, xcorr, ycorr, phi, eta);
      //      const double ref_radius = towergeom->get_radius();

      //      phi = 0;
      //      if (phi > M_PI) phi -= 2. * M_PI;  // convert to [-pi,pi]]

      prob = -1;
      chi2 = 0;
      ndf = 0;
      prob = pp->GetProb(chi2, ndf);
      //      printf("Prob/Chi2/NDF= %f %f %d Ecl=%f\n",prob,chi2,ndf,ecl);

      cluster = new RawClusterv1();
      cluster->set_energy(ecl);
      cluster->set_ecore(ecore);

      cluster->set_r(sqrt(xg * xg + yg * yg));
      cluster->set_phi(atan2(yg, xg));
      cluster->set_z(zg);

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
        iy = ich / NBINX;
        ix = ich % NBINX;
        // that code needs a closer look - here are the towers
        // with their energy added to the cluster object where
        // the id is the tower id
        // !!!!! Make sure twrkey is correctly extracted
	//        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(towers->getCalorimeterID(), ix + BINX0, iy + BINY0);
        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(towers->getCalorimeterID(), iy + BINY0, ix + BINX0); // Becuase in this part index1 is iy
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
