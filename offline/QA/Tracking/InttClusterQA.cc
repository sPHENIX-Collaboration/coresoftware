#include "InttClusterQA.h"

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
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
InttClusterQA::InttClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int InttClusterQA::InitRun(PHCompositeNode * /*unused*/)
{
  for (const auto &layer : {0, 1, 2, 3})
  {
    if (layer < 2)
    {
      m_layerLadderMap.insert(std::make_pair(layer, 12));
    }
    else
    {
      m_layerLadderMap.insert(std::make_pair(layer, 16));
    }
  }

  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttClusterQA::process_event(PHCompositeNode *topNode)
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

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::inttId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto ladderphiid = InttDefs::getLadderPhiId(hsk);
    auto sensor = InttDefs::getLadderZId(hsk);

    if (m_sensorInfo)
    {
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto cluskey = iter->first;
        auto *const cluster = iter->second;
        auto globalpos = tGeometry->getGlobalPosition(cluskey, cluster);
        auto phi = atan2(globalpos(1), globalpos(0));
        auto clayer = TrkrDefs::getLayer(cluskey);
        h_cluspersensor[(int) (layer) -3][(int) ladderphiid][(int) sensor]->Fill(cluster->getLocalY(), cluster->getLocalX());
        h_clusPhi_incl->Fill(phi);
        if (clayer == 3 || clayer == 4)
        {
          h_clusPhi_l34->Fill(phi);
          h_clusZ_clusPhi_l34->Fill(globalpos(2), phi);
        }
        else if (clayer == 5 || clayer == 6)
        {
          h_clusPhi_l56->Fill(phi);
          h_clusZ_clusPhi_l56->Fill(globalpos(2), phi);
        }

        m_totalClusters++;
        numclusters++;
      }
      m_nclustersPerSensor[((int) layer) - 3][(int) ladderphiid][(int) sensor] += numclusters;
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
        h_clusPhi_incl->Fill(phi);
        if (clayer == 3 || clayer == 4)
        {
          h_clusPhi_l34->Fill(phi);
          h_clusZ_clusPhi_l34->Fill(globalpos(2), phi);
        }
        else if (clayer == 5 || clayer == 6)
        {
          h_clusPhi_l56->Fill(phi);
          h_clusZ_clusPhi_l56->Fill(globalpos(2), phi);
        }
      }
    }
  }

  TrkrHitSetContainer::ConstRange hitsetrange = trkrHitSetContainer->getHitSets(TrkrDefs::TrkrId::inttId);

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    int sensor_hits = hitsetitr->second->size();
    float sensor_occupancy = (float) sensor_hits / (128. * 26.);
    h_occupancy->Fill(100. * sensor_occupancy);
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}
int InttClusterQA::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string InttClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void InttClusterQA::createHistos()
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  h_occupancy = new TH1F(std::format("{}sensorOccupancy", getHistoPrefix()).c_str(), "INTT Sensor Occupancy", 100, 0, 5);
  h_occupancy->GetXaxis()->SetTitle("Sensor Occupancy [%]");
  h_occupancy->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_occupancy);
  h_clusSize = new TH1F(std::format("{}clusterSize", getHistoPrefix()).c_str(), "INTT Cluster Size", 20, -0.5, 19.5);
  h_clusSize->GetXaxis()->SetTitle("Cluster Size");
  h_clusSize->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusSize);
  h_clusPhi_incl = new TH1F(std::format("{}clusterPhi_incl", getHistoPrefix()).c_str(), "INTT Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_incl->GetXaxis()->SetTitle("Cluster (inner+outer) #phi [rad]");
  h_clusPhi_incl->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_incl);
  h_clusPhi_l34 = new TH1F(std::format("{}clusterPhi_l34", getHistoPrefix()).c_str(), "INTT Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_l34->GetXaxis()->SetTitle("Cluster (inner) #phi [rad]");
  h_clusPhi_l34->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l34);
  h_clusPhi_l56 = new TH1F(std::format("{}clusterPhi_l56", getHistoPrefix()).c_str(), "INTT Cluster Phi", 320, -3.2, 3.2);
  h_clusPhi_l56->GetXaxis()->SetTitle("Cluster (outer) #phi [rad]");
  h_clusPhi_l56->GetYaxis()->SetTitle("Entries");
  hm->registerHisto(h_clusPhi_l56);
  h_clusZ_clusPhi_l34 = new TH2F(std::format("{}clusterZ_clusPhi_l34", getHistoPrefix()).c_str(), "INTT Cluster Z vs Cluster Phi", 55, cluszbin, 350, -3.5, 3.5);
  h_clusZ_clusPhi_l34->GetXaxis()->SetTitle("Cluster (inner) Z [cm]");
  h_clusZ_clusPhi_l34->GetYaxis()->SetTitle("Cluster (inner) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l34);
  h_clusZ_clusPhi_l56 = new TH2F(std::format("{}clusterZ_clusPhi_l56", getHistoPrefix()).c_str(), "INTT Cluster Z vs Cluster Phi", 55, cluszbin, 350, -3.5, 3.5);
  h_clusZ_clusPhi_l56->GetXaxis()->SetTitle("Cluster (outer) Z [cm]");
  h_clusZ_clusPhi_l56->GetYaxis()->SetTitle("Cluster (outer) #phi [rad]");
  hm->registerHisto(h_clusZ_clusPhi_l56);

  if (m_sensorInfo)
  {
    for (const auto &[layer, ladders] : m_layerLadderMap)
    {
      for (int ladder = 0; ladder < ladders; ladder++)
      {
        //! 4 sensor on each ladder
        for (int sensor = 0; sensor < 4; sensor++)
        {
          h_cluspersensor[layer][ladder][sensor] = new TH2F(std::format("{}ncluspersensor{}_{}_{}", getHistoPrefix(), layer, ladder, sensor).c_str(),
                                                            "INTT clusters per sensor", 100, -5, 5, 1000, -1, 1);
          h_cluspersensor[layer][ladder][sensor]->GetXaxis()->SetTitle("Local z [cm]");
          h_cluspersensor[layer][ladder][sensor]->GetYaxis()->SetTitle("Local rphi [cm]");
          hm->registerHisto(h_cluspersensor[layer][ladder][sensor]);
        }
      }
    }
  }

  return;
}
