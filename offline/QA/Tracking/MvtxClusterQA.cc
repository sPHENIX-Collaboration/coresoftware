
#include "MvtxClusterQA.h"

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>

#include <format>

//____________________________________________________________________________..
MvtxClusterQA::MvtxClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MvtxClusterQA::InitRun(PHCompositeNode *topNode)
{
  auto *geomContainer = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geomContainer)
  {
    std::cout << PHWHERE
              << " CYLINDERGEOM_MVTX  node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  for (const auto &layer : {0, 1, 2})
  {
    auto *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geomContainer->GetLayerGeom(layer));
    if (!layergeom)
    {
      std::cout << PHWHERE << "Did not get layergeom for layer "
                << layer << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    int nstaves = layergeom->get_N_staves();
    m_layerStaveMap.insert(std::make_pair(layer, nstaves));
  }

  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxClusterQA::process_event(PHCompositeNode *topNode)
{
  auto *clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "No cluster container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto *trkrHitSetContainer = findNode::getClass<TrkrHitSetContainerv1>(topNode, "TRKR_HITSET");
  if (!trkrHitSetContainer)
  {
    std::cout << PHWHERE << "No trkrhitset container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto *tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts geometry on node tree, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int numclusters = 0;
  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
  {
    auto range = clusterContainer->getClusters(hsk);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      numclusters++;
    }
  }

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
  {
    int numclusters_chip = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto stave = MvtxDefs::getStaveId(hsk);
    auto chip = MvtxDefs::getChipId(hsk);

    if (m_chipInfo)
    {
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto cluskey = iter->first;
        auto *const cluster = iter->second;
        auto globalpos = tGeometry->getGlobalPosition(cluskey, cluster);
        auto phi = atan2(globalpos(1), globalpos(0));
        auto clayer = TrkrDefs::getLayer(cluskey);
        h_clusperchip[(int) layer][(int) stave][(int) chip]->Fill(cluster->getLocalY(), cluster->getLocalX());
        h_clusSize->Fill(cluster->getSize());
        h_clusSize_nClus->Fill(numclusters, cluster->getSize());
        h_clusPhi_incl->Fill(phi);
        if (clayer == 0)
        {
          h_clusPhi_l0->Fill(phi);
          h_clusZ_clusPhi_l0->Fill(globalpos(2), phi);
        }
        else if (clayer == 1)
        {
          h_clusPhi_l1->Fill(phi);
          h_clusZ_clusPhi_l1->Fill(globalpos(2), phi);
        }
        else if (clayer == 2)
        {
          h_clusPhi_l2->Fill(phi);
          h_clusZ_clusPhi_l2->Fill(globalpos(2), phi);
        }
        m_totalClusters++;
        numclusters_chip++;
      }

      m_nclustersPerChip[(int) layer][(int) stave][(int) chip] += numclusters_chip;
    }
    else
    {
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto cluskey = iter->first;
        auto *const cluster = iter->second;
        auto globalpos = tGeometry->getGlobalPosition(cluskey, cluster);
        auto phi = atan2(globalpos(1), globalpos(0));
        auto clayer = TrkrDefs::getLayer(cluskey);
        h_clusSize->Fill(cluster->getSize());
        h_clusSize_nClus->Fill(numclusters, cluster->getSize());
        h_clusPhi_incl->Fill(phi);
        if (clayer == 0)
        {
          h_clusPhi_l0->Fill(phi);
          h_clusZ_clusPhi_l0->Fill(globalpos(2), phi);
        }
        else if (clayer == 1)
        {
          h_clusPhi_l1->Fill(phi);
          h_clusZ_clusPhi_l1->Fill(globalpos(2), phi);
        }
        else if (clayer == 2)
        {
          h_clusPhi_l2->Fill(phi);
          h_clusZ_clusPhi_l2->Fill(globalpos(2), phi);
        }
      }
    }
  }

  TrkrHitSetContainer::ConstRange hitsetrange = trkrHitSetContainer->getHitSets(TrkrDefs::TrkrId::mvtxId);

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    int chip_hits = hitsetitr->second->size();
    float chip_occupancy = (float) chip_hits / (512 * 1024);
    chip_occupancy = 100 * chip_occupancy;
    h_occupancy->Fill(chip_occupancy);
    int strobe = MvtxDefs::getStrobeId(hitsetitr->first);
    for (int i = 0; i < chip_hits; i++)
    {
      h_strobe->Fill(strobe);
    }
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}
int MvtxClusterQA::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..

std::string MvtxClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void MvtxClusterQA::createHistos()
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  h_occupancy = new TH1F(std::format("{}chipOccupancy", getHistoPrefix()).c_str(), "MVTX Chip Occupancy", 60, 0, 0.6);
  h_occupancy->GetXaxis()->SetTitle("Chip Occupancy [%]");
  h_occupancy->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_occupancy);
  h_clusSize = new TH1F(std::format("{}clusterSize", getHistoPrefix()).c_str(), "MVTX Cluster Size", 50, -0.5, 49.5);
  h_clusSize->GetXaxis()->SetTitle("Cluster Size");
  h_clusSize->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusSize);
  h_clusPhi_incl = new TH1F(std::format("{}clusterPhi_incl", getHistoPrefix()).c_str(), "MVTX Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_incl->GetXaxis()->SetTitle("Cluster (layer 0+1+2) #phi [rad]");
  h_clusPhi_incl->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_incl);
  h_clusPhi_l0 = new TH1F(std::format("{}clusterPhi_l0", getHistoPrefix()).c_str(), "MVTX Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_l0->GetXaxis()->SetTitle("Cluster (layer 0) #phi [rad]");
  h_clusPhi_l0->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l0);
  h_clusPhi_l1 = new TH1F(std::format("{}clusterPhi_l1", getHistoPrefix()).c_str(), "MVTX Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_l1->GetXaxis()->SetTitle("Cluster (layer 1) #phi [rad]");
  h_clusPhi_l1->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l1);
  h_clusPhi_l2 = new TH1F(std::format("{}clusterPhi_l2", getHistoPrefix()).c_str(), "MVTX Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_l2->GetXaxis()->SetTitle("Cluster (layer 2) #phi [rad]");
  h_clusPhi_l2->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l2);
  h_clusZ_clusPhi_l0 = new TH2F(std::format("{}clusterZ_clusPhi_l0", getHistoPrefix()).c_str(), "MVTX Cluster Z vs Phi", 300, -15, 15, 350, -3.5, 3.5);
  h_clusZ_clusPhi_l0->GetXaxis()->SetTitle("Cluster (layer 2) Z [cm]");
  h_clusZ_clusPhi_l0->GetYaxis()->SetTitle("Cluster (layer 0) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l0);
  h_clusZ_clusPhi_l1 = new TH2F(std::format("{}clusterZ_clusPhi_l1", getHistoPrefix()).c_str(), "MVTX Cluster Z vs Phi", 300, -15, 15, 350, -3.5, 3.5);
  h_clusZ_clusPhi_l1->GetXaxis()->SetTitle("Cluster (layer 2) Z [cm]");
  h_clusZ_clusPhi_l1->GetYaxis()->SetTitle("Cluster (layer 1) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l1);
  h_clusZ_clusPhi_l2 = new TH2F(std::format("{}clusterZ_clusPhi_l2", getHistoPrefix()).c_str(), "MVTX Cluster Z vs Phi", 300, -15, 15, 350, -3.5, 3.5);
  h_clusZ_clusPhi_l2->GetXaxis()->SetTitle("Cluster (layer 2) Z [cm]");
  h_clusZ_clusPhi_l2->GetYaxis()->SetTitle("Cluster (layer 2) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l2);
  h_strobe = new TH1I(std::format("{}strobeTiming", getHistoPrefix()).c_str(), "MVTX Strobe Timing per Hit", 32, -15, 16);
  h_strobe->GetXaxis()->SetTitle("Strobe BCO - GL1 BCO");
  h_strobe->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_strobe);
  h_clusSize_nClus = new TH2F(std::format("{}clusSize_nCLus", getHistoPrefix()).c_str(), "MVTX Cluster Size vs Number of Clusters", 800, -0.5, 799.5, 25, -0.5, 24.5);
  h_clusSize_nClus->GetXaxis()->SetTitle("Number of Clusters");
  h_clusSize_nClus->GetYaxis()->SetTitle("Cluster Size");
  hm->registerHisto(h_clusSize_nClus);

  if (m_chipInfo)
  {
    for (const auto &[layer, nstave] : m_layerStaveMap)
    {
      for (int stave = 0; stave < nstave; stave++)
      {
        //! 9 chips on each stave
        for (int chip = 0; chip < 9; chip++)
        {
          h_clusperchip[layer][stave][chip] = new TH2F(std::format("{}nclusperchip{}_{}_{}", getHistoPrefix(), layer, stave, chip).c_str(),
                                                       "MVTX clusters per chip", 2000, -2, 2, 2000, -1, 1);
          h_clusperchip[layer][stave][chip]->GetXaxis()->SetTitle("Local z [cm]");
          h_clusperchip[layer][stave][chip]->GetYaxis()->SetTitle("Local rphi [cm]");
          hm->registerHisto(h_clusperchip[layer][stave][chip]);
        }
      }
    }
  }

  return;
}
