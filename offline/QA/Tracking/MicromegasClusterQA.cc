#include "MicromegasClusterQA.h"

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>

#include <boost/format.hpp>

//____________________________________________________________________________..
MicromegasClusterQA::MicromegasClusterQA(const std::string &name)
  : SubsysReco(name)
{
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
    auto h = dynamic_cast<TH2 *>(hm->getHisto((boost::format("%sncluspertile%i_%i") % getHistoPrefix() % (((int) layer) - 55) % (int) tile).str()));

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
int MicromegasClusterQA::EndRun(const int /*runnumber*/)
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
 
  for (const auto &[layer, ntiles] : m_layerTileMap)
  {
    for (int tile = 0; tile < ntiles; tile++)
    {
      auto h = new TH2F((boost::format("%sncluspertile%i_%i") % getHistoPrefix() % layer % tile).str().c_str(),
                        "Micromegas clusters per tile", 2000, -30, 30, 2000, -20, 20);
      h->GetXaxis()->SetTitle("Local z [cm]");
      h->GetYaxis()->SetTitle("Local rphi [cm]");
      hm->registerHisto(h);

    }
  }

  return;
}
