#include "RawClusterBuilderTemplate.h"

#include "BEmcCluster.h"
#include "BEmcRec.h"
#include "BEmcRecCEMC.h"

#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
// To isolate the issue with the vertex
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterv2.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>


#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

RawClusterBuilderTemplate::RawClusterBuilderTemplate(const std::string &name)
  : SubsysReco(name)
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
  else
  {
    std::cout << "Warning from RawClusterBuilderTemplate::Detector(): no detector specific class "
              << Name() << " defined for detector " << detector
              << ". Default BEmcRec will be used" << std::endl;
    bemc = new BEmcRec();
  }

  // Set vertex
  float vertex[3] = {0, 0, 0};
  bemc->SetVertex(vertex);
  // Set threshold
  bemc->SetTowerThreshold(_min_tower_e);
  bemc->SetPeakThreshold(_min_peak_e);
  bemc->SetProbNoiseParam(fProbNoiseParam);
}

void RawClusterBuilderTemplate::set_UseCorrectPosition(const bool useCorrectPosition)
{
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::set_UseCorrectPosition()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return;
  }

  bemc->set_UseCorrectPosition(useCorrectPosition);
}

void RawClusterBuilderTemplate::set_UseCorrectShowerDepth(const bool useCorrectShowerDepth) {
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::set_UseCorrectShowerDepth()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return;
  }

  bemc->set_UseCorrectShowerDepth(useCorrectShowerDepth);
}

void RawClusterBuilderTemplate::set_UseDetailedGeometry(const bool useDetailedGeometry)
{
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::set_UseDetailedGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return;
  }

  m_UseDetailedGeometry = useDetailedGeometry;
  bemc->set_UseDetailedGeometry(m_UseDetailedGeometry);
}

void RawClusterBuilderTemplate::LoadProfile(const std::string &fname)
{
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::LoadProfile()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return;
  }
  std::string url = CDBInterface::instance()->getUrl("EMCPROFILE", fname);
  bemc->LoadProfile(url);
}

void RawClusterBuilderTemplate::SetCylindricalGeometry()
{
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::SetCylindricalGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return;
  }

  bemc->SetCylindricalGeometry();
}

void RawClusterBuilderTemplate::SetPlanarGeometry()
{
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::SetPlanarGeometry()(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return;
  }

  bemc->SetPlanarGeometry();
}

int RawClusterBuilderTemplate::InitRun(PHCompositeNode *topNode)
{
  if (bemc == nullptr)
  {
    std::cerr << "Error in RawClusterBuilderTemplate::InitRun(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Ensure that the detailed geometry is available if the user requests it.
  // Otherwise, use the default geometry 
  if (m_UseDetailedGeometry && detector != "CEMC")
  {
    m_UseDetailedGeometry = false;
    bemc->set_UseDetailedGeometry(false);
    std::cout << "Warning in RawClusterBuilderTemplate::InitRun()(): No alternative detailed geometry defined for detector " << detector << ". m_UseDetailedGeometry automatically set to false." << std::endl;
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

  // geometry case
  if (m_TowerGeomNodeName.empty())
  {
    m_TowerGeomNodeName = "TOWERGEOM_" + detector;
    if (m_UseDetailedGeometry)
    {
      if (detector == "CEMC")
      {
        m_TowerGeomNodeName = m_TowerGeomNodeName + "_DETAILED";
      }
      else
      {
        std::cout << "RawClusterBuilderTemplate::InitRun - Detailed geometry not implemented for detector " << detector << ". The former geometry is used instead" << std::endl;
      }
    }
  }
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << m_TowerGeomNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int Calo_ID = towergeom->get_calorimeter_id();
  // std::cout << std::endl << std::endl << std::endl << "Calorimeter ID: " << Calo_ID << std::endl << std::endl << std::endl;

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
    ixmin = std::min(ixmin, ix);
    ixmax = std::max(ixmax, ix);
    iymin = std::min(iymin, iy);
    iymax = std::max(iymax, iy);
    ngeom++;
  }
  if (Verbosity() > 1)
  {
    std::cout << "Info from RawClusterBuilderTemplate::InitRun(): Init geometry for "
              << detector << ": N of geom towers: " << ngeom << "; ix = "
              << ixmin << "-" << ixmax << ", iy = "
              << iymin << "-" << iymax << std::endl;
  }
  if (ixmax < ixmin || iymax < iymin)
  {
    std::cout << "Error in RawClusterBuilderTemplate::InitRun(): wrong geometry data for detector "
              << detector << std::endl;
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

    // reconstruction change
    bemc->SetTowerGeometry(ix, iy, *towerg);
    bemc->SetCalotype(Calo_ID);
    if (Calo_ID == RawTowerDefs::EEMC ||
        Calo_ID == RawTowerDefs::EEMC_crystal ||
        Calo_ID == RawTowerDefs::EEMC_glass)
    {
      bemc->SetScinSize(towerg->get_size_z());
    }
  }


  // geometry case
  // With the former geometry, the tower dimensions are approximated from
  // the consecutive tower center postions in the "CompleteTowerGeometry" method.
  // With the detailed geometry, the information is already given in the 
  // RawTowerGeom node so no further step is required.
  if (!(m_UseDetailedGeometry && detector == "CEMC"))
  {
    if (!bemc->CompleteTowerGeometry())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // geometry case
  if (bPrintGeom)
  {
    std::string fname;
    if (m_UseDetailedGeometry && detector == "CEMC")
    {
      fname = "geom_" + detector + "_detailed.txt";
      bemc->PrintTowerGeometryDetailed(fname);
    }
    else
    {
      fname = "geom_" + detector + ".txt";
      bemc->PrintTowerGeometry(fname);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTemplate::PrintCylGeom(RawTowerGeomContainer *towergeom, const std::string &fname) const
{
  std::ofstream outfile(fname);
  if (!outfile.is_open())
  {
    std::cout << "Error in BEmcRec::RawClusterBuilderTemplate::PrintCylGeom(): Failed to open file "
              << fname << std::endl;
    return;
  }
  outfile << NBINX << " " << NBINY << std::endl;
  for (int ip = 0; ip < NBINX; ip++)
  {
    outfile << ip << " " << towergeom->get_phicenter(ip) << std::endl;
  }
  for (int ip = 0; ip < NBINY; ip++)
  {
    outfile << ip << " " << towergeom->get_etacenter(ip) << std::endl;
  }
  outfile.close();
}

bool RawClusterBuilderTemplate::Cell2Abs(RawTowerGeomContainer * /*towergeom*/, float /*phiC*/, float /*etaC*/, float &phi, float &eta)
{
  phi = eta = 0;
  return false;
}

int RawClusterBuilderTemplate::process_event(PHCompositeNode *topNode)
{
  if (bemc == nullptr)
  {
    std::cout << "Error in RawClusterBuilderTemplate::process_event(): detector is not defined; use RawClusterBuilderTemplate::Detector() to define it" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawTowerContainer *towers = nullptr;
  if (m_UseTowerInfo < 1)
  {
    std::string towernodename = "TOWER_CALIB_" + detector;

    // Grab the towers
    towers = findNode::getClass<RawTowerContainer>(topNode, towernodename);
    if (!towers)
    {
      std::cout << PHWHERE << ": Could not find node " << towernodename << std::endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  }

  m_TowerGeomNodeName = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << m_TowerGeomNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  TowerInfoContainer *calib_towerinfos = nullptr;
  if (m_UseTowerInfo > 0)
  {
    std::string towerinfoNodename = "TOWERINFO_CALIB_" + detector;
    if (!m_inputnodename.empty())
    {
      towerinfoNodename = m_inputnodename;
    }

    calib_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerinfoNodename);
    if (!calib_towerinfos)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << " " << towerinfoNodename << " Node missing, doing bail out!"
                << std::endl;

      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  }

  // Get vertex
  float vx = 0;
  float vy = 0;
  float vz = 0;
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (vertexmap && m_UseAltZVertex == 0)  // default
  {
    GlobalVertex* vtx = vertexmap->begin()->second;
    if (vtx)
    {
      auto typeStartIter = vtx->find_vertexes(m_vertex_type);
      auto typeEndIter = vtx->end_vertexes();
      for (auto iter = typeStartIter; iter != typeEndIter; ++iter)
      {
        const auto& [type, vertexVec] = *iter;
        if (type != m_vertex_type)
        {
          continue;
        }
        for (const auto* vertex : vertexVec)
        {
          if (!vertex)
          {
            continue;
          }
          vx = vertex->get_x();
          vy = vertex->get_y();
          vz = vertex->get_z();
        }
      }
    }
  }

  MbdVertexMap *mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");

  if (mbdmap && m_UseAltZVertex == 1)
  {
    MbdVertex *bvertex = nullptr;
    for (MbdVertexMap::ConstIter mbditer = mbdmap->begin();
         mbditer != mbdmap->end();
         ++mbditer)
    {
      bvertex = mbditer->second;
    }
    if (bvertex)
    {
      vz = bvertex->get_z();
    }
  }
  
  if (m_UseAltZVertex == 3)
  {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {
      PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetVtxRange();
      for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter) 
      {
         PHG4VtxPoint *vtx_tr = iter->second;
         if ( vtx_tr->get_id() == 1 )
         { 
           vz = vtx_tr->get_z();
           vy = vtx_tr->get_y();
           vx = vtx_tr->get_x();
         }
      }
    }
    else {
      std::cout << "RawClusterBuilderTemplate: Error requiring truth vertex but non was found. Exiting" << std::endl;  
      return Fun4AllReturnCodes::ABORTEVENT;
    } 
  }

  // Set vertex
  float vertex[3] = {vx, vy, vz};
  bemc->SetVertex(vertex);
  // Set threshold
  bemc->SetTowerThreshold(_min_tower_e);
  bemc->SetPeakThreshold(_min_peak_e);

  bemc->SetProbNoiseParam(fProbNoiseParam);
  bemc->SetProfileProb(bProfProb);

  _clusters->Reset();  // make sure cluster container is empty before filling it with new clusters

  // Define vector of towers in EmcModule format to input into BEmc
  EmcModule vhit;
  std::vector<EmcModule> HitList;
  HitList.erase(HitList.begin(), HitList.end());
  int ich;

  if (m_UseTowerInfo < 1)
  {
    // make the list of towers above threshold
    RawTowerContainer::ConstRange begin_end = towers->getTowers();
    RawTowerContainer::ConstIterator itr = begin_end.first;

    for (; itr != begin_end.second; ++itr)
    {
      RawTower *tower = itr->second;
      // debug change
      //      std::cout << "  Tower e = " << tower->get_energy()
      //           << " (" << _min_tower_e << ")" << std::endl;
      if (IsAcceptableTower(tower))
      {
        // debug change
        // std::cout << "(" << tower->get_column() << "," << tower->get_row()
        //      << ")  (" << tower->get_binphi() << "," << tower->get_bineta()
        //      << ")" << std::endl;
        //	  ix = tower->get_column();
        RawTowerDefs::keytype towerid = tower->get_id();
        int ix = RawTowerDefs::decode_index2(towerid);  // index2 is phi in CYL
        int iy = RawTowerDefs::decode_index1(towerid);  // index1 is eta in CYL
        ix -= BINX0;
        iy -= BINY0;
        // debug change
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
  }
  else if (m_UseTowerInfo)
  {
    // make the list of towers above threshold
    // debug change
    // TowerInfoContainer::ConstRange begin_end = calib_towerinfos->getTowers();
    // TowerInfoContainer::ConstIterator rtiter;
    unsigned int nchannels = calib_towerinfos->size();
    // for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)

    //float total_energy = 0;
    for (unsigned int channel = 0; channel < nchannels; channel++)
    {
      TowerInfo *tower_info = calib_towerinfos->get_tower_at_channel(channel);
      
      // debug change
      //      std::cout << "  Tower e = " << tower->get_energy()
      //           << " (" << _min_tower_e << ")" << std::endl;
      if (IsAcceptableTower(tower_info))
      {
        unsigned int towerkey = calib_towerinfos->encode_key(channel);
        int ieta = calib_towerinfos->getTowerEtaBin(towerkey);
        int iphi = calib_towerinfos->getTowerPhiBin(towerkey);

        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);

        int ix = RawTowerDefs::decode_index2(key);  // index2 is phi in CYL
        int iy = RawTowerDefs::decode_index1(key);  // index1 is eta in CYL

        ix -= BINX0;
        iy -= BINY0;

        if (ix >= 0 && ix < NBINX && iy >= 0 && iy < NBINY)
        {
          ich = iy * NBINX + ix;
          // add key field to vhit
          vhit.ich = ich;
          vhit.amp = tower_info->get_energy() * fEnergyNorm;  // !!! Global Calibration
          vhit.tof = tower_info->get_time();

          // debug change
          /*total_energy += vhit.amp;
          
          if (vhit.amp > 1e-8) {
            std::cout << "hit: (" << iy << ", " << ix << ", " << vhit.amp << ", " << tower_info->get_energy() << ", " << vhit.tof << ")\n";
          }
          */

          HitList.push_back(vhit);
        }
      }
    }
    // debug change
    //std::cout << "Total hit energy = " << total_energy << std::endl;
  }

  bemc->SetModules(&HitList);

  // Find clusters (as a set of towers with common edge)
  int ncl = bemc->FindClusters();
  if (ncl < 0)
  {
    std::cout << "!!! Error in BEmcRec::FindClusters(): numbers of cluster "
              << ncl << " ?" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get pointer to clusters
  std::vector<EmcCluster> *ClusterList = bemc->GetClusters();
  std::vector<EmcCluster>::iterator pc;

  std::vector<EmcCluster>::iterator pp;
  float ecl;
  float ecore;
  float xcg;
  float ycg;
  float xx;
  float xy;
  float yy;
  //  float xcorr, ycorr;
  EmcModule hmax;
  RawCluster *cluster;

  std::vector<EmcCluster> PList;
  std::vector<EmcModule> Peaks;

  float prob;
  float chi2;
  int ndf;
  float xg;
  float yg;
  float zg;

  std::vector<EmcModule>::iterator ph;
  std::vector<EmcModule> hlist;

  // ncl = 0;
  for (pc = ClusterList->begin(); pc != ClusterList->end(); ++pc)
  {
    //    ecl = pc->GetTotalEnergy();
    //    pc->GetMoments( &xcg, &ycg, &xx, &xy, &yy );

    int npk = pc->GetSubClusters(PList, Peaks,m_subclustersplitting);
    if (npk < 0)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    // debug change
    //    std::cout << "  iCl = " << ncl << " (" << npk << "): E ="
    //         << ecl << "  x = " << xcg << "  y = " << ycg << std::endl;

    for (pp = PList.begin(); pp != PList.end(); ++pp)
    {
      // Cluster energy
      ecl = pp->GetTotalEnergy();
      if (ecl < m_min_cluster_e)
      {
        continue;
      }
      ecore = pp->GetECoreCorrected();
      // 3x3 energy around center of gravity
      // e9 = pp->GetE9();
      // Ecore (basically near 2x2 energy around center of gravity)
      // ecore = pp->GetECore();
      // Center of Gravity etc.
      pp->GetMoments(xcg, ycg, xx, xy, yy);

      if (m_UseAltZVertex == 2)
      {
        xg = -99999;  // signal to force zvtx = 0
        pp->GetGlobalPos(xg, yg, zg);
      }
      else
      {
        xg = 0;  // usual mode, uses regular zvtx
        pp->GetGlobalPos(xg, yg, zg);
      }

      // Tower with max energy
      hmax = pp->GetMaxTower();

      // debug change
      //      phi = (xcg-float(NPHI)/2.+0.5)/float(NPHI)*2.*M_PI;
      //      eta = (ycg-float(NETA)/2.+0.5)/float(NETA)*2.2; // -1.1<eta<1.1;

      //      Cell2Abs(towergeom,xcg,ycg,phi,eta);

      //      pp->GetCorrPos(&xcorr, &ycorr);
      //      Cell2Abs(towergeom, xcorr, ycorr, phi, eta);
      //      const double ref_radius = towergeom->get_radius();

      //      phi = 0;
      //      if (phi > M_PI) phi -= 2. * M_PI;  // convert to [-pi,pi]]

      //      prob = -1;
      chi2 = 0;
      ndf = 0;
      prob = pp->GetProb(chi2, ndf);
      //      std::cout << "Prob/Chi2/NDF = " << prob << " " << chi2
      //           << " " << ndf << " Ecl = " << ecl << std::endl;

      cluster = m_writeClusterV2
                    ? static_cast<RawCluster*>(new RawClusterv2())
                    : static_cast<RawCluster*>(new RawClusterv1());
      cluster->set_energy(ecl);
      cluster->set_ecore(ecore);
      cluster->set_r(std::sqrt(xg * xg + yg * yg));
      cluster->set_phi(std::atan2(yg, xg));
      cluster->set_z(zg);

      cluster->set_prob(prob);
      if (ndf > 0)
      {
        cluster->set_chi2(chi2 / ndf);
      }
      else
      {
        cluster->set_chi2(0);
      }
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
        const RawTowerDefs::CalorimeterId Calo_ID = towergeom->get_calorimeter_id();
        RawTowerDefs::keytype twrkey = RawTowerDefs::encode_towerid(Calo_ID, iy + BINY0, ix + BINX0);  // Becuase in this part index1 is iy
        //	std::cout << iphi << " " << ieta << ": "
        //           << twrkey << " e = " << (*ph).amp) << std::endl;
        cluster->addTower(twrkey, (*ph).amp / fEnergyNorm);
        ++ph;
      }

      // stamp tower CoG (raw & corrected) only when writing v2
      if (m_writeClusterV2)
      {
        float xcorr = xcg, ycorr = ycg;
        bemc->CorrectPosition(ecl, xcg, ycg, xcorr, ycorr);
        if (auto* c2 = dynamic_cast<RawClusterv2*>(cluster))
        {
          c2->set_tower_cog(xcg, ycg, xcorr, ycorr);
        }
      }

      auto it_v2same = _clusters->AddCluster(cluster);
      cluster->set_id(it_v2same->first);
      // ncl++;

      //      std::cout << "    ipk = " << ipk << ": E = " << ecl << "  E9 = "
      //           << e9 << "  x = " << xcg << "  y = " << ycg
      //           << "  MaxTower: (" << hmax.ich%NPHI << ","
      //           << hmax.ich/NPHI << ") e = " << hmax.amp << std::endl;
    }
  }

  if (chkenergyconservation && towers && _clusters)
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
  else if (chkenergyconservation)
  {
    std::cout << "RawClusterBuilderTemplate : energy conservation check asked for but tower or cluster container is NULL" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTemplate::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the cEMC node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  // Get the _det_name subnode
  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", detector));

  // Check that it is there
  if (!cemcNode)
  {
    cemcNode = new PHCompositeNode(detector);
    dstNode->addNode(cemcNode);
  }
  ClusterNodeName = "CLUSTER_" + detector;
  if (!m_outputnodename.empty())
  {
    ClusterNodeName = m_outputnodename;
  }
  else if (m_UseTowerInfo)
  {
    ClusterNodeName = "CLUSTERINFO_" + detector;
  }
  _clusters = findNode::getClass<RawClusterContainer>(dstNode, ClusterNodeName);
  if (!_clusters)
  {
    _clusters = new RawClusterContainer();
  }

  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName, "PHObject");
  cemcNode->addNode(clusterNode);
}

bool RawClusterBuilderTemplate::IsAcceptableTower(TowerInfo *tower) const
{
  if (tower->get_energy() < _min_tower_e)
  {
    return false;
  }

  if (m_do_tower_selection)
  {
    if (!tower->get_isGood())
    {
      return false;
    }

  }
  return true;
}

bool RawClusterBuilderTemplate::IsAcceptableTower(RawTower *tower) const
{
  if (tower->get_energy() < _min_tower_e)
  {
    return false;
  }
  return true;
}
