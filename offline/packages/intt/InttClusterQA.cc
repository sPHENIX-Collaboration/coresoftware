#include "InttClusterQA.h"
#include "CylinderGeomIntt.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
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
InttClusterQA::InttClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int InttClusterQA::InitRun(PHCompositeNode * /*unused*/)
{
  for (auto &layer : {0, 1, 2, 3})
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

  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::inttId))
  {
    int numclusters = 0;
    auto range = clusterContainer->getClusters(hsk);
    auto layer = TrkrDefs::getLayer(hsk);
    auto ladderphiid = InttDefs::getLadderPhiId(hsk);
    auto sensor = InttDefs::getLadderZId(hsk);
    auto h = dynamic_cast<TH2 *>(hm->getHisto( boost::str(boost::format("%sncluspersensor%i_%i_%i") %getHistoPrefix() %(((int) layer) - 3) %((int) ladderphiid) %((int) sensor)).c_str()));

    for (auto iter = range.first; iter != range.second; ++iter)
    {
      // const auto cluskey = iter->first;
      const auto cluster = iter->second;
      h->Fill(cluster->getLocalY(), cluster->getLocalX());
      m_totalClusters++;
      numclusters++;
    }
    m_nclustersPerSensor[((int) layer) - 3][(int) ladderphiid][(int) sensor] += numclusters;
  }

  m_event++;

  return Fun4AllReturnCodes::EVENT_OK;
}
int InttClusterQA::EndRun(const int runnumber)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2 *h_totalclusters = dynamic_cast<TH2 *>(hm->getHisto( boost::str(boost::format("%snclusperrun") %getHistoPrefix()).c_str()));
  // NOLINTNEXTLINE(bugprone-integer-division)
  h_totalclusters->Fill(runnumber, m_totalClusters / m_event);

  for (const auto &[layer, ladders] : m_layerLadderMap)
  {
    for (int ladder = 0; ladder < ladders; ladder++)
    {
      for (int sensor = 0; sensor < 4; sensor++)
      {
        TH2 *h = dynamic_cast<TH2 *>(hm->getHisto( boost::str(boost::format("%sncluspersensorperrun%i_%i_%i") %getHistoPrefix() %layer %ladder %sensor).c_str()));
        if (h)
        {
          // NOLINTNEXTLINE(bugprone-integer-division)
          h->Fill(runnumber, m_nclustersPerSensor[layer][ladder][sensor] / m_event);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string InttClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}
void InttClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto h = new TH2F( boost::str(boost::format("%snclusperrun") %getHistoPrefix()).c_str(),
                      "INTT Clusters per event per run number", m_runbins, m_beginRun, m_endRun, 1000, 0, 1000);
    h->GetXaxis()->SetTitle("Run number");
    h->GetYaxis()->SetTitle("Clusters per event");
    hm->registerHisto(h);
  }

  for (const auto &[layer, ladders] : m_layerLadderMap)
  {
    for (int ladder = 0; ladder < ladders; ladder++)
    {
      //! 4 sensor on each ladder
      for (int sensor = 0; sensor < 4; sensor++)
      {
        auto h = new TH2F( boost::str(boost::format("%sncluspersensor%i_%i_%i") %getHistoPrefix()
				      %layer %ladder %sensor).c_str(),
                          "INTT clusters per sensor", 100, -5, 5, 1000, -1, 1);
        h->GetXaxis()->SetTitle("Local z [cm]");
        h->GetYaxis()->SetTitle("Local rphi [cm]");
        hm->registerHisto(h);

	auto h2 = new TH2F( boost::str(boost::format("%sncluspersensorperrun%i_%i_%i") %getHistoPrefix() %layer %ladder %sensor).c_str(),
                           "INTT clusters per event per sensor per run", m_runbins, m_beginRun, m_endRun, 100, 0, 100);
        h2->GetXaxis()->SetTitle("Run number");
        h2->GetYaxis()->SetTitle("Clusters per event");
        hm->registerHisto(h2);
      }
    }
  }

  return;
}
