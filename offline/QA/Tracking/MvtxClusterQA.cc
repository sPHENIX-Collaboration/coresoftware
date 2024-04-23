
#include "MvtxClusterQA.h"
#include <mvtx/CylinderGeom_Mvtx.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <qautils/QAHistManagerDef.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/MvtxDefs.h>
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
MvtxClusterQA::MvtxClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MvtxClusterQA::Init(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxClusterQA::InitRun(PHCompositeNode *topNode)
{
  auto geomContainer = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geomContainer)
  {
    std::cout << PHWHERE
              << " CYLINDERGEOM_MVTX  node not found on node tree"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  for (auto &layer : {0, 1, 2})
  {
    auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geomContainer->GetLayerGeom(layer));
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

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto stave = MvtxDefs::getStaveId(hsk);
    auto chip = MvtxDefs::getChipId(hsk);
    auto h = dynamic_cast<TH2 *>(hm->getHisto(Form("%snclusperchip%i_%i_%i", getHistoPrefix().c_str(), (int) layer, (int) stave, (int) chip)));

    for (auto iter = range.first; iter != range.second; ++iter)
    {
      // const auto cluskey = iter->first;
      const auto cluster = iter->second;
      h->Fill(cluster->getLocalY(), cluster->getLocalX());
      m_totalClusters++;
      numclusters++;
    }
    m_nclustersPerChip[(int) layer][(int) stave][(int) chip] += numclusters;
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}
int MvtxClusterQA::EndRun(const int runnumber)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2 *h_totalclusters = dynamic_cast<TH2 *>(hm->getHisto(Form("%snclusperrun", getHistoPrefix().c_str())));
  h_totalclusters->Fill(runnumber, m_totalClusters / m_event);

  for (const auto &[layer, staves] : m_layerStaveMap)
  {
    for (int stave = 0; stave < staves; stave++)
    {
      for (int chip = 0; chip < 9; chip++)
      {
        TH2 *h = dynamic_cast<TH2 *>(hm->getHisto(Form("%snclusperchipperrun%i_%i_%i", getHistoPrefix().c_str(), layer, stave, chip)));
        if (h)
        {
          h->Fill(runnumber, m_nclustersPerChip[layer][stave][chip] / m_event);
        }
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int MvtxClusterQA::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string MvtxClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}
void MvtxClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto h = new TH2F(Form("%snclusperrun", getHistoPrefix().c_str()),
                      "MVTX Clusters per event per run number", m_runbins, m_beginRun, m_endRun, 1000, 0, 1000);
    h->GetXaxis()->SetTitle("Run number");
    h->GetYaxis()->SetTitle("Clusters per event");
    hm->registerHisto(h);
  }

  for (const auto &[layer, nstave] : m_layerStaveMap)
  {
    for (int stave = 0; stave < nstave; stave++)
    {
      //! 9 chips on each stave
      for (int chip = 0; chip < 9; chip++)
      {
        auto h = new TH2F(Form("%snclusperchip%i_%i_%i", getHistoPrefix().c_str(),
                               layer, stave, chip),
                          "MVTX clusters per chip", 2000, -2, 2, 2000, -1, 1);
        h->GetXaxis()->SetTitle("Local z [cm]");
        h->GetYaxis()->SetTitle("Local rphi [cm]");
        hm->registerHisto(h);

        auto h2 = new TH2F(Form("%snclusperchipperrun%i_%i_%i", getHistoPrefix().c_str(), layer, stave, chip),
                           "MVTX clusters per event per chip per run", m_runbins, m_beginRun, m_endRun, 100, 0, 100);
        h2->GetXaxis()->SetTitle("Run number");
        h2->GetYaxis()->SetTitle("Clusters per event");
        hm->registerHisto(h2);
      }
    }
  }

  return;
}
