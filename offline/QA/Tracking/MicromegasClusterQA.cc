
#include "MicromegasClusterQA.h"
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <qautils/QAHistManagerDef.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <TH1F.h>
#include <TH2F.h>

//____________________________________________________________________________..
MicromegasClusterQA::MicromegasClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
MicromegasClusterQA::~MicromegasClusterQA()
{
}

//____________________________________________________________________________..
int MicromegasClusterQA::Init(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasClusterQA::InitRun(PHCompositeNode *topNode)
{
  auto geomContainer = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!geomContainer)
  {
    std::cout << PHWHERE
              << " CYLINDERGEOM_MICROMEGAS_FULL  node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  int layer = 0;
  const auto range = geomContainer->get_begin_end();

  for (auto iter = range.first; iter != range.second; ++iter)
  {
    auto layergeom = static_cast<CylinderGeomMicromegas *>(iter->second);
    int ntiles = layergeom->get_tiles_count();
    m_layerTileMap.insert(std::make_pair(layer, ntiles));
    layer++;
  }

  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasClusterQA::process_event(PHCompositeNode *topNode)
{
  auto clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "No cluster container, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts geometry on node tree, bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::micromegasId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto tile = MicromegasDefs::getTileId(hsk);
    auto h = dynamic_cast<TH2 *>(hm->getHisto(Form("%sncluspertile%i_%i", getHistoPrefix().c_str(), ((int) layer) - 55, (int) tile)));

    for (auto iter = range.first; iter != range.second; ++iter)
    {
      // const auto cluskey = iter->first;
      const auto cluster = iter->second;
      h->Fill(cluster->getLocalY(), cluster->getLocalX());
      m_totalClusters++;
      numclusters++;
    }
    m_nclustersPerTile[((int) layer) - 55][(int) tile] += numclusters;
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}
int MicromegasClusterQA::EndRun(const int runnumber)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2 *h_totalclusters = dynamic_cast<TH2 *>(hm->getHisto(Form("%snclusperrun", getHistoPrefix().c_str())));
  h_totalclusters->Fill(runnumber, m_totalClusters / m_event);

  for (const auto &[layer, tiles] : m_layerTileMap)
  {
    for (int tile = 0; tile < tiles; tile++)
    {
      TH2 *h = dynamic_cast<TH2 *>(hm->getHisto(Form("%sncluspertileperrun%i_%i", getHistoPrefix().c_str(), layer, tile)));
      if (h)
      {
        h->Fill(runnumber, m_nclustersPerTile[layer][tile] / m_event);
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int MicromegasClusterQA::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string MicromegasClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}
void MicromegasClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto h = new TH2F(Form("%snclusperrun", getHistoPrefix().c_str()),
                      "Micromegas Clusters per event per run number", m_runbins, m_beginRun, m_endRun, 1000, 0, 1000);
    h->GetXaxis()->SetTitle("Run number");
    h->GetYaxis()->SetTitle("Clusters per event");
    hm->registerHisto(h);
  }

  for (const auto &[layer, ntiles] : m_layerTileMap)
  {
    for (int tile = 0; tile < ntiles; tile++)
    {
      auto h = new TH2F(Form("%sncluspertile%i_%i", getHistoPrefix().c_str(),
                             layer, tile),
                        "Micromegas clusters per tile", 2000, -30, 30, 2000, -20, 20);
      h->GetXaxis()->SetTitle("Local z [cm]");
      h->GetYaxis()->SetTitle("Local rphi [cm]");
      hm->registerHisto(h);

      auto h2 = new TH2F(Form("%sncluspertileperrun%i_%i", getHistoPrefix().c_str(), layer, tile),
                         "Micromegas clusters per event per tile per run", m_runbins, m_beginRun, m_endRun, 100, 0, 20);
      h2->GetXaxis()->SetTitle("Run number");
      h2->GetYaxis()->SetTitle("Clusters per event");
      hm->registerHisto(h2);
    }
  }

  return;
}
