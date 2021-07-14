#include "RawClusterBuilderTemplate.h"

#include "BEmcCluster.h"
#include "BEmcRec.h"
#include "BEmcRecCEMC.h"
#include "BEmcRecEEMC.h"
#include "BEmcRecFEMC.h"

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
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

#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

RawClusterBuilderTemplate::RawClusterBuilderTemplate(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
  , bemc(nullptr)
  , fEnergyNorm(1.)
  , _min_tower_e(0.020)
  , chkenergyconservation(0)
  , detector("NONE")
  , BINX0(0)
  , NBINX(0)
  , BINY0(0)
  , NBINY(0)
  , bPrintGeom(false)
  , bProfProb(false)
{
}

RawClusterBuilderTemplate::~RawClusterBuilderTemplate()
{
  // one can delete null pointers
  delete bemc;
}

void RawClusterBuilderTemplate::Detector(const std::string &d)
{
  detector = d;

  // Create proper BEmcRec object

  if (detector == "CEMC")
  {
    bemc = new BEmcRecCEMC();
  }
  else if (detector == "FEMC")
  {
    bemc = new BEmcRecFEMC();
  }
  else if (detector == "EEMC")
  {
    bemc = new BEmcRecEEMC();
  }
  else if (detector == "EEMC_crystal")
  {
    bemc = new BEmcRecEEMC();
  }
  else if (detector == "EEMC_glass")
  {
    bemc = new BEmcRecEEMC();
  }
  else
  {
    cout << "Warning from RawClusterBuilderTemplate::Detector(): no detector specific class "
         << Name() << " defined for detector " << detector
         << ". Default BEmcRec will be used" << endl;
    bemc = new BEmcRec();
  }

  // Set vertex
  float vertex[3] = {0, 0, 0};
  bemc->SetVertex(vertex);
  // Set threshold
  bemc->SetTowerThreshold(0);
}

void RawClusterBuilderTemplate::LoadProfile(const string &fname)
{
  //  _emcprof = new BEmcProfile(fname);
  bemc->LoadProfile(fname);
}

void RawClusterBuilderTemplate::SetCylindricalGeometry()
{
  if (bemc == nullptr)
  {
    cout << "Error in RawClusterBuilderTemplate::SetCylindricalGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << endl;
    return;
  }

  bemc->SetCylindricalGeometry();
}

void RawClusterBuilderTemplate::SetPlanarGeometry()
{
  if (bemc == nullptr)
  {
    cout << "Error in RawClusterBuilderTemplate::SetPlanarGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << endl;
    return;
  }

  bemc->SetPlanarGeometry();
}

int RawClusterBuilderTemplate::InitRun(PHCompositeNode *topNode)
{
  if (bemc == nullptr)
  {
    cout << "Error in RawClusterBuilderTemplate::InitRun(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << endl;
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
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!towergeom)
  {
    cout << PHWHERE << ": Could not find node " << towergeomnodename << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int Calo_ID = towergeom->get_calorimeter_id();
  // cout << endl << endl << endl << "Calorimeter ID: " << Calo_ID << endl << endl << endl;

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
    int ix = RawTowerDefs::decode_index2(towerid);  // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(towerid);  // index1 is eta in CYL
    if (ixmin > ix) ixmin = ix;
    if (ixmax < ix) ixmax = ix;
    if (iymin > iy) iymin = iy;
    if (iymax < iy) iymax = iy;
    ngeom++;
  }
  if (Verbosity() > 1)
  {
    cout << "Info from RawClusterBuilderTemplate::InitRun(): Init geometry for "
         << detector << ": N of geom towers: " << ngeom << "; ix = "
         << ixmin << "-" << ixmax << ", iy = "
         << iymin << "-" << iymax << endl;
  }
  if (ixmax < ixmin || iymax < iymin)
  {
    cout << "Error in RawClusterBuilderTemplate::InitRun(): wrong geometry data for detector "
         << detector << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  BINX0 = ixmin;
  NBINX = ixmax - ixmin + 1;
  BINY0 = iymin;
  NBINY = iymax - iymin + 1;

  bemc->SetDim(NBINX, NBINY);

  itr_geom = begin_end_geom.first;
  for (; itr_geom != begin_end_geom.second; ++itr_geom)
  {
    RawTowerGeom *towerg = itr_geom->second;
    RawTowerDefs::keytype towerid = towerg->get_id();
    //    int itype = towerg->get_tower_type();
    //    if( itype==2 ) { // PbSc
    int ix = RawTowerDefs::decode_index2(towerid);  // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(towerid);  // index1 is eta in CYL
    ix -= BINX0;
    iy -= BINY0;

    bemc->SetTowerGeometry(ix, iy, towerg->get_center_x(), towerg->get_center_y(), towerg->get_center_z());
    bemc->SetCalotype(Calo_ID);
    if (Calo_ID == RawTowerDefs::EEMC ||
        Calo_ID == RawTowerDefs::EEMC_crystal ||
        Calo_ID == RawTowerDefs::EEMC_glass)
    {
      bemc->SetScinSize(towerg->get_size_z());
    }
  }

  if (!bemc->CompleteTowerGeometry()) return Fun4AllReturnCodes::ABORTEVENT;

  if (bPrintGeom)
  {
    string fname = "geom_" + detector + ".txt";
    //    bemc->PrintTowerGeometry("geom.txt");
    bemc->PrintTowerGeometry(fname);
    //    PrintCylGeom(towergeom,"phieta.txt");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTemplate::PrintCylGeom(RawTowerGeomContainer *towergeom, const string &fname)
{
  ofstream outfile(fname);
  if (!outfile.is_open())
  {
    cout << "Error in BEmcRec::RawClusterBuilderTemplate::PrintCylGeom(): Failed to open file "
         << fname << endl;
    return;
  }
  outfile << NBINX << " " << NBINY << endl;
  for (int ip = 0; ip < NBINX; ip++)
  {
    outfile << ip << " " << towergeom->get_phicenter(ip) << endl;
  }
  for (int ip = 0; ip < NBINY; ip++)
  {
    outfile << ip << " " << towergeom->get_etacenter(ip) << endl;
  }
  outfile.close();
}

bool RawClusterBuilderTemplate::Cell2Abs(RawTowerGeomContainer *towergeom, float phiC, float etaC, float &phi, float &eta)
{
  phi = eta = 0;
  return false;
}

int RawClusterBuilderTemplate::process_event(PHCompositeNode *topNode)
{
  if (bemc == nullptr)
  {
    cout << "Error in RawClusterBuilderTemplate::process_event(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  string towernodename = "TOWER_CALIB_" + detector;
  // Grab the towers
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernodename);
  if (!towers)
  {
    std::cout << PHWHERE << ": Could not find node " << towernodename << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  string towergeomnodename = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!towergeom)
  {
    cout << PHWHERE << ": Could not find node " << towergeomnodename << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get vertex
  float vx = 0;
  float vy = 0;
  float vz = 0;
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (vertexmap)
  {
    if (!vertexmap->empty())
    {
      GlobalVertex *vertex = (vertexmap->begin()->second);
      vx = vertex->get_x();
      vy = vertex->get_y();
      vz = vertex->get_z();
    }
  }

  // Set vertex
  float vertex[3] = {vx, vy, vz};
  bemc->SetVertex(vertex);
  // Set threshold
  bemc->SetTowerThreshold(_min_tower_e);

  bemc->SetProfileProb(bProfProb);

  // _clusters->Reset(); // !!! Not sure if it is necessarry to do it - ask Chris

  // make the list of towers above threshold
  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator itr = begin_end.first;

  // Define vector of towers in EmcModule format to input into BEmc
  EmcModule vhit;
  std::vector<EmcModule> HitList;
  HitList.erase(HitList.begin(), HitList.end());
  int ich;

  for (; itr != begin_end.second; ++itr)
  {
    RawTower *tower = itr->second;
    //      cout << "  Tower e = " << tower->get_energy()
    //           << " (" << _min_tower_e << ")" << endl;
    if (tower->get_energy() > _min_tower_e)
    {
      // cout << "(" << tower->get_column() << "," << tower->get_row()
      //      << ")  (" << tower->get_binphi() << "," << tower->get_bineta()
      //      << ")" << endl;
      //	  ix = tower->get_column();
      RawTowerDefs::keytype towerid = tower->get_id();
      int ix = RawTowerDefs::decode_index2(towerid);  // index2 is phi in CYL
      int iy = RawTowerDefs::decode_index1(towerid);  // index1 is eta in CYL
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
    cout << "!!! Error in BEmcRec::FindClusters(): numbers of cluster "
         << ncl << " ?" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get pointer to clusters
  std::vector<EmcCluster> *ClusterList = bemc->GetClusters();
  std::vector<EmcCluster>::iterator pc;

  std::vector<EmcCluster>::iterator pp;
  float ecl, ecore, xcg, ycg, xx, xy, yy;
  //  float xcorr, ycorr;
  EmcModule hmax;
  RawCluster *cluster;

  std::vector<EmcCluster> PList;
  std::vector<EmcModule> Peaks;
  std::vector<EmcCluster> *pPList = &PList;
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

    int npk = pc->GetSubClusters(pPList, pPeaks);
    if (npk < 0) return Fun4AllReturnCodes::ABORTEVENT;

    //    cout << "  iCl = " << ncl << " (" << npk << "): E ="
    //         << ecl << "  x = " << xcg << "  y = " << ycg << endl;

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
      pp->GetMoments(xcg, ycg, xx, xy, yy);
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
      //      cout << "Prob/Chi2/NDF = " << prob << " " << chi2
      //           << " " << ndf << " Ecl = " << ecl << endl;

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
        int iy = ich / NBINX;
        int ix = ich % NBINX;
        // that code needs a closer look - here are the towers
        // with their energy added to the cluster object where
        // the id is the tower id
        // !!!!! Make sure twrkey is correctly extracted
        //        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(towers->getCalorimeterID(), ix + BINX0, iy + BINY0);
        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(towers->getCalorimeterID(), iy + BINY0, ix + BINX0);  // Becuase in this part index1 is iy
        //	cout << iphi << " " << ieta << ": "
        //           << twrkey << " e = " << (*ph).amp) << endl;
        cluster->addTower(twrkey, (*ph).amp / fEnergyNorm);
        ++ph;
      }

      _clusters->AddCluster(cluster);
      // ncl++;

      //      cout << "    ipk = " << ipk << ": E = " << ecl << "  E9 = "
      //           << e9 << "  x = " << xcg << "  y = " << ycg
      //           << "  MaxTower: (" << hmax.ich%NPHI << ","
      //           << hmax.ich/NPHI << ") e = " << hmax.amp << endl;
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
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName, "PHObject");
  cemcNode->addNode(clusterNode);
}
