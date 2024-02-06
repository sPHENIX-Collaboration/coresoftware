#include "RawClusterPositionCorrection.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>  // for decode_index1, decode_in...
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <cdbobjects/CDBHistos.h>  // for CDBHistos
#include <cdbobjects/CDBTTree.h>   // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>

#include <algorithm>  // for max

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>  // for pair

RawClusterPositionCorrection::RawClusterPositionCorrection(const std::string &name)
  : SubsysReco(std::string("RawClusterPositionCorrection_") + name)
  , _det_name(name)
  , bins_eta(384)
  , bins_phi(64)
  , iEvent(0)
{
}
RawClusterPositionCorrection::~RawClusterPositionCorrection()
{
  delete cdbHisto;
  delete cdbttree;
}

int RawClusterPositionCorrection::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);

  // access the cdb and get cdbtree
  std::string m_calibName = "cemc_PDC_NorthSouth_8x8_23instru";
  std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
  if (!calibdir.empty())
  {
    cdbttree = new CDBTTree(calibdir);
  }
  else
  {
    std::cout << std::endl
              << "did not find CDB tree" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  h2NorthSector = new TH2F("h2NorthSector", "Cluster; towerid #eta; towerid #phi", bins_eta, 47.5, 95.5, bins_phi, -0.5, 7.5);
  h2SouthSector = new TH2F("h2SouthSector", "Cluster; towerid #eta; towerid #phi", bins_eta, -0.5, 47.5, bins_phi, -0.5, 7.5);

  /// north
  std::string m_fieldname = "cemc_PDC_NorthSector_8x8_clusE";
  std::string m_fieldname_ecore = "cemc_PDC_NorthSector_8x8_clusEcore";
  float calib_constant = 0;

  // Read in the calibration factors and store in the array
  for (int i = 0; i < bins_phi; ++i)
  {
    std::vector<float> dumvec;
    std::vector<float> dumvec2;
    for (int j = 0; j < bins_eta; ++j)
    {
      int key = i * bins_eta + j;
      calib_constant = cdbttree->GetFloatValue(key, m_fieldname);
      dumvec.push_back(calib_constant);
      calib_constant = cdbttree->GetFloatValue(key, m_fieldname_ecore);
      dumvec2.push_back(calib_constant);
    }
    calib_constants_north.push_back(dumvec);
    calib_constants_north_ecore.push_back(dumvec2);
  }
  /// south
  m_fieldname = "cemc_PDC_SouthSector_8x8_clusE";
  m_fieldname_ecore = "cemc_PDC_SouthSector_8x8_clusEcore";

  // Read in the calibration factors and store in the array
  for (int i = 0; i < bins_phi; ++i)
  {
    std::vector<float> dumvec;
    std::vector<float> dumvec2;
    for (int j = 0; j < bins_eta; ++j)
    {
      int key = i * bins_eta + j;
      calib_constant = cdbttree->GetFloatValue(key, m_fieldname);
      dumvec.push_back(calib_constant);
      calib_constant = cdbttree->GetFloatValue(key, m_fieldname_ecore);
      dumvec2.push_back(calib_constant);
    }
    calib_constants_south.push_back(dumvec);
    calib_constants_south_ecore.push_back(dumvec2);
  }

  // Load PDC final stage correction
  calibdir = CDBInterface::instance()->getUrl("cemc_PDC_ResidualCorr");
  if (!calibdir.empty())
  {
    cdbHisto = new CDBHistos(calibdir.c_str());
    cdbHisto->LoadCalibrations();
    //pdcCorrFlat = cdbHisto->getHisto("h1_res_p");
    pdcCorrFlat = cdbHisto->getHisto("h_res_E_eta");
  }
  else
  {
    std::cout << std::endl
              << "did not find CDB histo" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterPositionCorrection::process_event(PHCompositeNode *topNode)
{
  // make sure new cluster container is empty 
  _recalib_clusters->Reset();

  if (Verbosity() >= Fun4AllBase::VERBOSITY_SOME)
  {
    if (iEvent % 100 == 0) std::cout << "Progress: " << iEvent << std::endl;
    ++iEvent;
  }

  std::string rawClusNodeName = "CLUSTER_" + _det_name;
  if (m_UseTowerInfo)
  {
    rawClusNodeName = "CLUSTERINFO_" + _det_name;
  }

  RawClusterContainer *rawclusters = findNode::getClass<RawClusterContainer>(topNode, rawClusNodeName.c_str());
  if (!rawclusters)
  {
    std::cout << "No " << _det_name << " Cluster Container found while in RawClusterPositionCorrection, can't proceed!!!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawTowerContainer *_towers = nullptr;
  TowerInfoContainer *_towerinfos = nullptr;

  if (!m_UseTowerInfo)
  {
    _towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_" + _det_name);
    if (!_towers)
    {
      std::cout << "No calibrated " << _det_name << " tower info found while in RawClusterPositionCorrection, can't proceed!!!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else
  {
    std::string towerinfoNodename = "TOWERINFO_CALIB_" + _det_name;

    _towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerinfoNodename);
    if (!_towerinfos)
    {
      std::cerr << Name() << "::" << _det_name << "::" << __PRETTY_FUNCTION__
                << " " << towerinfoNodename << " Node missing, doing bail out!"
                << std::endl;

      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  }

  std::string towergeomnodename = "TOWERGEOM_" + _det_name;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnodename << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  const int nphibin = towergeom->get_phibins();

  // loop over the clusters
  RawClusterContainer::ConstRange begin_end = rawclusters->getClusters();
  RawClusterContainer::ConstIterator iter;

  for (iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    //    RawClusterDefs::keytype key = iter->first;
    RawCluster *cluster = iter->second;

    float clus_energy = cluster->get_energy();

    std::vector<float> toweretas;
    std::vector<float> towerphis;
    std::vector<float> towerenergies;

    // loop over the towers in the cluster
    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;
    for (toweriter = towers.first;
         toweriter != towers.second;
         ++toweriter)
    {
      if (!m_UseTowerInfo)
      {
        RawTower *tower = _towers->getTower(toweriter->first);

        int towereta = tower->get_bineta();
        int towerphi = tower->get_binphi();
        double towerenergy = tower->get_energy();

        // put the etabin, phibin, and energy into the corresponding vectors
        toweretas.push_back(towereta);
        towerphis.push_back(towerphi);
        towerenergies.push_back(towerenergy);
      }
      else
      {
        // std::cout << "clus " << iter->first << " tow key " << toweriter->first << " decoded to " << towerindex << std::endl;

        int iphi = RawTowerDefs::decode_index2(toweriter->first);  // index2 is phi in CYL
        int ieta = RawTowerDefs::decode_index1(toweriter->first);  // index1 is eta in CYL
        unsigned int towerkey = iphi + (ieta << 16U);

        assert(_towerinfos);

        unsigned int towerindex = _towerinfos->decode_key(towerkey);

        TowerInfo *towinfo = _towerinfos->get_tower_at_channel(towerindex);

        double towerenergy = towinfo->get_energy();

        // put the eta, phi, energy into corresponding vectors
        toweretas.push_back(ieta);
        towerphis.push_back(iphi);
        towerenergies.push_back(towerenergy);
      }
    }

    int ntowers = toweretas.size();
    //    std::cout << "jf " <<  ntowers << std::endl;
    assert(ntowers >= 1);

    // loop over the towers to determine the energy
    // weighted eta and phi position of the cluster

    float etamult = 0;
    float etasum = 0;
    float phimult = 0;
    float phisum = 0;

    for (int j = 0; j < ntowers; j++)
    {
      float energymult = towerenergies.at(j) * toweretas.at(j);
      etamult += energymult;
      etasum += towerenergies.at(j);

      int phibin = towerphis.at(j);

      if (phibin - towerphis.at(0) < -nphibin / 2.0)
      {
        phibin += nphibin;
      }
      else if (phibin - towerphis.at(0) > +nphibin / 2.0)
      {
        phibin -= nphibin;
      }
      assert(std::abs(phibin - towerphis.at(0)) <= nphibin / 2.0);

      energymult = towerenergies.at(j) * phibin;
      phimult += energymult;
      phisum += towerenergies.at(j);
    }

    float avgphi = phimult / phisum;
    if (isnan(avgphi)) continue;

    float avgeta = etamult / etasum;

    if (avgphi < 0)
    {
      avgphi += nphibin;
    }

    avgphi = fmod(avgphi, nphibin);

    if (avgphi >= 255.5) avgphi -= bins_phi;

    avgphi = fmod(avgphi + 0.5, 8) - 0.5;  // wrapping [-0.5, 255.5] to [-0.5, 7.5]
    int etabin = -99;
    int phibin = -99;

    // check if the cluster is in the north or south sector
    if (avgeta < 47.5)
    {
      etabin = h2SouthSector->GetXaxis()->FindBin(avgeta) - 1;
    }
    else
    {
      etabin = h2NorthSector->GetXaxis()->FindBin(avgeta) - 1;
    }
    phibin = h2NorthSector->GetYaxis()->FindBin(avgphi) - 1;  // can use either h2NorthSector or h2SouthSector since both have the same phi binning

    if ((phibin < 0 || etabin < 0) && Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
    {
      std::cout << "couldn't recalibrate cluster, something went wrong??" << std::endl;
    }

    float eclus_recalib_val = 1;
    float ecore_recalib_val = 1;
    if (phibin > -1 && etabin > -1)
    {
      if (avgeta < 47.5)
      {
        eclus_recalib_val = calib_constants_south[phibin][etabin];
        ecore_recalib_val = calib_constants_south_ecore[phibin][etabin];
      }
      else
      {
        eclus_recalib_val = calib_constants_north[phibin][etabin];
        ecore_recalib_val = calib_constants_north_ecore[phibin][etabin];
      }
    }
    RawCluster *recalibcluster = dynamic_cast<RawCluster *>(cluster->CloneMe());

    assert(recalibcluster);
    //    if (m_UseTowerInfo)
    recalibcluster->set_energy(clus_energy / eclus_recalib_val);
    recalibcluster->set_ecore(cluster->get_ecore() / ecore_recalib_val);

    CLHEP::Hep3Vector vertex(0,0,0);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recalibcluster, vertex);
    float clusEta = E_vec_cluster.pseudoRapidity();

    if (cluster->get_ecore()    >= pdcCorrFlat->GetXaxis()->GetXmin() 
        && cluster->get_ecore() <  pdcCorrFlat->GetXaxis()->GetXmax()
        && clusEta              >= pdcCorrFlat->GetYaxis()->GetXmin()
        && clusEta              <  pdcCorrFlat->GetYaxis()->GetXmax()
       )
    {

      int ecoreBin = pdcCorrFlat->GetXaxis()->FindBin(recalibcluster->get_ecore());
      int etaBin = pdcCorrFlat -> GetYaxis() -> FindBin(clusEta);
      float pdcCalib =  pdcCorrFlat -> GetBinContent(ecoreBin, etaBin);
      //float pdcCalib = pdcCorrFlat->GetBinContent(ecoreBin);
      if (pdcCalib < 0.1) pdcCalib = 1;

      recalibcluster->set_ecore(recalibcluster->get_ecore() / pdcCalib);
    }

    _recalib_clusters->AddCluster(recalibcluster);

    if (Verbosity() >= Fun4AllBase::VERBOSITY_EVEN_MORE && clus_energy > 1)
    {
      std::cout << "Input eclus cluster energy: " << clus_energy << std::endl;
      std::cout << "Recalib value: " << eclus_recalib_val << std::endl;
      std::cout << "phibin: " << phibin << ", etabin: " << etabin << std::endl;
      std::cout << "Recalibrated eclus cluster energy: "
                << clus_energy / eclus_recalib_val << std::endl;
      std::cout << "Input ecore cluster energy: "
                << cluster->get_ecore() << std::endl;
      std::cout << "Recalib value: " << ecore_recalib_val << std::endl;
      std::cout << "Recalibrated ecore cluster energy: "
                << cluster->get_ecore() / ecore_recalib_val << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


void RawClusterPositionCorrection::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Get the DST Node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // Check that it is there
  if (!dstNode)
  {
    std::cout << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in RawClusterPositionCorrection::CreateNodeTree");
  }

  // Get the _det_name subnode
  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _det_name));

  // Check that it is there
  if (!cemcNode)
  {
    cemcNode = new PHCompositeNode(_det_name);
    dstNode->addNode(cemcNode);
  }

  // Check to see if the cluster recalib node is on the nodetree
  std::string ClusterCorrNodeName = "CLUSTER_POS_COR_" + _det_name;
  if (m_UseTowerInfo)
  {
    ClusterCorrNodeName = "CLUSTERINFO_POS_COR_" + _det_name;
  }
  _recalib_clusters = findNode::getClass<RawClusterContainer>(topNode, ClusterCorrNodeName);
  if (_recalib_clusters)
  {
    _recalib_clusters->Clear();
  }
  else
  {
    _recalib_clusters = new RawClusterContainer();
    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_recalib_clusters, ClusterCorrNodeName, "PHObject");
    cemcNode->addNode(clusterNode);
  }
}
int RawClusterPositionCorrection::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
