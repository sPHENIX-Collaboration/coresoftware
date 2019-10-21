#include "RawClusterBuilderTemplateFEMC.h"

#include "BEmcCluster.h"
#include "BEmcProfile.h"
#include "BEmcRecFEMC.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <cstdio>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

RawClusterBuilderTemplateFEMC::~RawClusterBuilderTemplateFEMC()
{
  delete bemc;
  delete _emcprof;
}

RawClusterBuilderTemplateFEMC::RawClusterBuilderTemplateFEMC(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , _emcprof(nullptr)
  , _min_tower_e(0.020)
  , chkenergyconservation(0)
  , detector("NONE")
{
  BINX0 = 0;
  NBINX = 0;
  BINY0 = 0;
  NBINY = 0;
  Zcenter = 0;

  fEnergyNorm = 1.;
  bemc = new BEmcRecFEMC();
  //
  // Some initial values for clustering
  //
  bemc->SetPlaneGeometry();
  // Configuration: number of towers in Phi and Eta, and tower dimension (here still in angle units and eta units)
  bemc->SetGeometry(64, 64, 1.0, 1.0);
  // Define vertex ... not used now
  float vertex[3] = {0, 0, 0};
  bemc->SetVertex(vertex);
  // Define threshold ... not used now
  bemc->SetTowerThreshold(0);
}

void RawClusterBuilderTemplateFEMC::LoadProfile(const char *fname)
{
  _emcprof = new BEmcProfile(fname);
}

int RawClusterBuilderTemplateFEMC::InitRun(PHCompositeNode *topNode)
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

  int ngeom = 0;
  int ixmin = 999999;
  int ixmax = -999999;
  int iymin = 999999;
  int iymax = -999999;
  float sz = 0;
  RawTowerGeomContainer::ConstRange begin_end_geom = towergeom->get_tower_geometries();
  RawTowerGeomContainer::ConstIterator itr_geom = begin_end_geom.first;
  for (; itr_geom != begin_end_geom.second; ++itr_geom)
  {
    RawTowerGeom *towerg = itr_geom->second;
    RawTowerDefs::keytype towerid = towerg->get_id();
    //    int itype = towerg->get_tower_type();
    //    if( itype==2 ) { // PbSc
    int ix = RawTowerDefs::decode_index1(towerid);
    int iy = RawTowerDefs::decode_index2(towerid);
    if (ixmin > ix) ixmin = ix;
    if (ixmax < ix) ixmax = ix;
    if (iymin > iy) iymin = iy;
    if (iymax < iy) iymax = iy;
    sz += towerg->get_center_z();
    ngeom++;
    //    }
  }
  Zcenter = 305.;
  if (ngeom > 0) Zcenter = sz / ngeom;

  printf("************* Init FEMC: N of geom towers: %d; ix=%d-%d iy=%d-%d Zcenter=%f\n", ngeom, ixmin, ixmax, iymin, iymax, Zcenter);

  if (ixmax < ixmin || iymax < iymin) return Fun4AllReturnCodes::ABORTEVENT;

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
    int ix = RawTowerDefs::decode_index1(towerid);
    int iy = RawTowerDefs::decode_index2(towerid);
    ix -= BINX0;
    iy -= BINY0;
    bemc->SetTowerGeometry(ix, iy, towerg->get_center_x(), towerg->get_center_y(), towerg->get_center_z());
    //    }
  }

  //  bemc->PrintTowerGeometry("geom_femc.txt");

  return Fun4AllReturnCodes::EVENT_OK;
}

bool RawClusterBuilderTemplateFEMC::Cell2Abs(RawTowerGeomContainer *towergeom, float xC, float yC, float &xA, float &yA)
{
  int ix = xC + 0.5;  // tower #
  if (ix < 0 || ix >= NBINX)
  {
    printf("RawClusterBuilderTemplateFEMC::Cell2Abs: wrong input x: %d\n", ix);
    return false;
  }

  int iy = yC + 0.5;  // tower #
  if (iy < 0 || iy >= NBINY)
  {
    printf("RawClusterBuilderTemplateFEMC::Cell2Abs: wrong input y: %d\n", iy);
    return false;
  }

  ix += BINX0;
  iy += BINY0;

  RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(towergeom->get_calorimeter_id(), ix, iy);
  RawTowerGeom *geom0 = towergeom->get_tower_geometry(key);

  // Next tower in x
  key = RawTowerDefs::encode_towerid(towergeom->get_calorimeter_id(), ix + 1, iy);
  RawTowerGeom *geomx = towergeom->get_tower_geometry(key);
  if (geomx == nullptr)
  {
    key = RawTowerDefs::encode_towerid(towergeom->get_calorimeter_id(), ix - 1, iy);
    geomx = towergeom->get_tower_geometry(key);
    if (geomx == nullptr)
    {
      printf("RawClusterBuilderTemplateFEMC::Cell2Abs: error in geometry extraction: %d %d\n", ix, iy);
      return false;
    }
  }

  // Next tower in y
  key = RawTowerDefs::encode_towerid(towergeom->get_calorimeter_id(), ix, iy + 1);
  RawTowerGeom *geomy = towergeom->get_tower_geometry(key);
  if (geomy == nullptr)
  {
    key = RawTowerDefs::encode_towerid(towergeom->get_calorimeter_id(), ix, iy - 1);
    geomy = towergeom->get_tower_geometry(key);
    if (geomy == nullptr)
    {
      printf("RawClusterBuilderTemplateFEMC::Cell2Abs: error in geometry extraction: %d %d\n", ix, iy);
      return false;
    }
  }

  //  float dx = geom->get_size_x();
  //  float dy = geom->get_size_y();
  float dx = fabs(geom0->get_center_x() - geomx->get_center_x());
  float dy = fabs(geom0->get_center_y() - geomy->get_center_y());
  xA = geom0->get_center_x() + (xC - ix + BINX0) * dx;
  yA = geom0->get_center_y() + (yC - iy + BINY0) * dy;

  return true;
}

int RawClusterBuilderTemplateFEMC::process_event(PHCompositeNode *topNode)
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
      ix = tower->get_bineta() - BINX0;  // eta: index1
      iy = tower->get_binphi() - BINY0;  // phi: index2
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

  float xout, yout, zout;
  float theta;
  float prob;
  float chi2 = NAN;
  int ndf;

  vector<EmcModule>::iterator ph;
  vector<EmcModule> hlist;

  //ncl = 0;
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

      //      pp->GetCorrPos(&xcorr,&ycorr);
      /* 
      xcorr = xcg;
      ycorr = ycg;
      Cell2Abs(towergeom,xcorr,ycorr,xout,yout);
      */
      //      const double ref_radius = towergeom->get_radius();

      pp->GetGlobalPos(xout, yout, zout);

      hlist = pp->GetHitList();
      float zVert = 0;  // !!!!! In future it should take actual zVert
      theta = atan(sqrt(xout * xout + yout * yout) / fabs(zout - zVert));

      //      prob = pp->GetProb(chi2,ndf);
      prob = -1;
      ndf = 0;
      if (_emcprof != nullptr) prob = _emcprof->GetProb(&hlist, NBINX, ecl, theta);

      cluster = new RawClusterv1();
      cluster->set_energy(ecl);
      cluster->set_ecore(ecore);

      cluster->set_r(sqrt(xout * xout + yout * yout));
      cluster->set_phi(atan2(yout, xout));
      cluster->set_z(zout);

      cluster->set_prob(prob);
      if (ndf > 0)
        cluster->set_chi2(chi2 / ndf);
      else
        cluster->set_chi2(0);

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
        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(towers->getCalorimeterID(), ix + BINX0, iy + BINY0);
        //	printf("%d %d: %d e=%f\n",iphi,ieta,twrkey,(*ph).amp);
        cluster->addTower(twrkey, (*ph).amp / fEnergyNorm);
        ++ph;
      }

      _clusters->AddCluster(cluster);
      //      ncl++;

      //      printf("    ipk=%d: E=%f x=%f (%f)  y=%f (%f)  MaxTower: (%d,%d) e=%f\n",ncl-1,ecl,xcorr,xout,ycorr,yout,hmax.ich%NBINX,hmax.ich/NBINX,hmax.amp);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderTemplateFEMC::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTemplateFEMC::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
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
