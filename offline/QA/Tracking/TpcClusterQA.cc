#include "TpcClusterQA.h"

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
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
TpcClusterQA::TpcClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcClusterQA::InitRun(PHCompositeNode *topNode)
{
  // find tpc geometry
  auto geomContainer =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geomContainer)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // TPC has 3 regions, inner, mid and outer
  std::vector<int> region_layer_low = {7, 23, 39};
  std::vector<int> region_layer_high = {22, 38, 54};

  // make a layer to region multimap
  const auto range = geomContainer->get_begin_end();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    m_layers.insert(iter->first);

    for (int region = 0; region < 3; ++region)
    {
      if (iter->first >= region_layer_low[region] && iter->first <= region_layer_high[region])
      {
        m_layerRegionMap.insert(std::make_pair(iter->first, region));
      }
    }
  }

  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcClusterQA::process_event(PHCompositeNode *topNode)
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

  TH2 *h_totalclusters = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "stotal_clusters")));
  TH2 *h_clusterssector = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "ncluspersector")));

  struct HistoList
  {
    TH1 *crphisize = nullptr;
    TH1 *czsize = nullptr;
    TH1 *crphierr = nullptr;
    TH1 *czerr = nullptr;
    TH1 *cedge = nullptr;
    TH1 *coverlap = nullptr;
  };

  int hitsetkeynum = 0;
  using HistoMap = std::map<int, HistoList>;
  HistoMap histos;

  for (auto &region : {0, 1, 2})
  {
    HistoList hist;
    hist.crphisize = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sphisize_%i") % getHistoPrefix() % region).str()));
    hist.czsize = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%szsize_%i") % getHistoPrefix() % region).str()));
    hist.crphierr = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%srphi_error_%i") % getHistoPrefix() % region).str()));
    hist.czerr = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sz_error_%i") % getHistoPrefix() % region).str()));
    hist.cedge = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusedge_%i") % getHistoPrefix() % region).str()));
    hist.coverlap = dynamic_cast<TH1 *>(hm->getHisto((boost::format("%sclusoverlap_%i") % getHistoPrefix() % region).str()));

    histos.insert(std::make_pair(region, hist));
  }
  auto fill = [](TH1 *h, float val)
  { if (h) { h->Fill(val); 
} };

  float nclusperevent[24] = {0};
  for (auto &hsk : clusterContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    int numclusters = 0;

    auto range = clusterContainer->getClusters(hsk);
    int sector = TpcDefs::getSectorId(hsk);
    int side = TpcDefs::getSide(hsk);
    if (side > 0)
    {
      sector += 12;
    }
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      const auto cluskey = iter->first;
      const auto cluster = iter->second;
      const auto it = m_layerRegionMap.find(TrkrDefs::getLayer(cluskey));
      int region = it->second;
      const auto hiter = histos.find(region);
      if (hiter == histos.end())
      {
        continue;
      }

      fill(hiter->second.crphisize, cluster->getPhiSize());
      fill(hiter->second.czsize, cluster->getZSize());
      fill(hiter->second.crphierr, cluster->getRPhiError());
      fill(hiter->second.czerr, cluster->getZError());
      fill(hiter->second.cedge, cluster->getEdge());
      fill(hiter->second.coverlap, cluster->getOverlap());

      numclusters++;
    }

    nclusperevent[sector] += numclusters;
    h_totalclusters->Fill(hitsetkeynum, numclusters);
    m_totalClusters += numclusters;
    hitsetkeynum++;
  }
  for (int i = 0; i < 24; i++)
  {
    h_clusterssector->Fill(i, nclusperevent[i]);
    m_clustersPerSector[i] += nclusperevent[i];
  }
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterQA::EndRun(const int runnumber)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2 *h_totalclusters = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nclusperrun")));
  h_totalclusters->Fill(runnumber, (float) m_totalClusters / m_event);

  for (int i = 0; i < 24; i++)
  {
    TH2 *h = dynamic_cast<TH2 *>(hm->getHisto((boost::format("%snclusperrun_sector%i") % getHistoPrefix() % i).str()));
    h->Fill(runnumber, (float) m_clustersPerSector[i] / m_event);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string TpcClusterQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}
void TpcClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nclusperrun").c_str(),
                      "TPC Clusters per event per run number", m_runbins, m_beginRun, m_endRun, 1000, 0, 1000);
    h->GetXaxis()->SetTitle("Run number");
    h->GetYaxis()->SetTitle("Clusters per event");
    hm->registerHisto(h);
  }
  {
    for (int i = 0; i < 24; i++)
    {
      auto h = new TH2F((boost::format("%snclusperrun_sector%i") % getHistoPrefix() % i).str().c_str(),
                        (boost::format("TPC Clusters per event per run number sector %i") % i).str().c_str(), m_runbins, m_beginRun, m_endRun, 1000, 0, 1000);
      h->GetXaxis()->SetTitle("Run number");
      h->GetYaxis()->SetTitle((boost::format("Clusters per event in Sector %i") % i).str().c_str());
      hm->registerHisto(h);
    }
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ncluspersector").c_str(),
                      "TPC Clusters per event per sector", 24, 0, 24, 1000, 0, 1000);
    h->GetXaxis()->SetTitle("Sector number");
    h->GetYaxis()->SetTitle("Clusters per event");
    hm->registerHisto(h);
  }
  for (auto &region : {0, 1, 2})
  {
    {
      auto h = new TH1F((boost::format("%sphisize_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC cluster #phi size region_%i") % region).str().c_str(), 10, 0, 10);
      h->GetXaxis()->SetTitle("Cluster #phi_{size}");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%szsize_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC cluster z size region_%i") % region).str().c_str(), 10, 0, 10);
      h->GetXaxis()->SetTitle("Cluster z_{size}");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%srphi_error_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC r#Delta#phi error region_%i") % region).str().c_str(), 100, 0, 0.075);
      h->GetXaxis()->SetTitle("r#Delta#phi error [cm]");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sz_error_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC z error region_%i") % region).str().c_str(), 100, 0, 0.18);
      h->GetXaxis()->SetTitle("z error [cm]");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusedge_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC hits on edge_%i") % region).str().c_str(), 30, 0, 30);
      h->GetXaxis()->SetTitle("Cluster edge");
      hm->registerHisto(h);
    }
    {
      auto h = new TH1F((boost::format("%sclusoverlap_%i") % getHistoPrefix() % region).str().c_str(),
                        (boost::format("TPC clus overlap_%i") % region).str().c_str(), 30, 0, 30);
      h->GetXaxis()->SetTitle("Cluster overlap");
      hm->registerHisto(h);
    }
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "stotal_clusters").c_str(),
                      "TPC clusters per hitsetkey", 1152, 0, 1152, 10000, 0, 10000);
    h->GetXaxis()->SetTitle("Hitsetkey number");
    h->GetYaxis()->SetTitle("Number of clusters");
    hm->registerHisto(h);
  }

  return;
}
