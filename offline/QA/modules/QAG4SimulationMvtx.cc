#include "QAG4SimulationMvtx.h"
#include "QAG4Util.h"
#include "QAHistManagerDef.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>  // for getTrkrId, getHit...
#include <trackbase/TrkrHitTruthAssoc.h>

#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TAxis.h>  // for TAxis
#include <TH1.h>
#include <TString.h>  // for Form

#include <cassert>
#include <cmath>     // for atan2
#include <iostream>  // for operator<<, basic...
#include <iterator>  // for distance
#include <map>       // for map
#include <utility>   // for pair, make_pair

//________________________________________________________________________
QAG4SimulationMvtx::QAG4SimulationMvtx(const std::string& name)
  : SubsysReco(name)
{
}

//________________________________________________________________________
int QAG4SimulationMvtx::InitRun(PHCompositeNode* topNode)
{
  // prevent multiple creations of histograms
  if (m_initialized)
    return Fun4AllReturnCodes::EVENT_OK;
  else
    m_initialized = true;

  // find mvtx geometry
  auto geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERGEOM_MVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get layers from mvtx geometry
  const auto range = geom_container->get_begin_end();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    m_layers.insert(iter->first);
  }

  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create histograms
  for (const auto& layer : m_layers)
  {
    if (Verbosity()) std::cout << PHWHERE << " adding layer " << layer << std::endl;
    {
      // rphi residuals (cluster - truth)
      auto h = new TH1F(Form("%sdrphi_%i", get_histo_prefix().c_str(), layer), Form("MVTX r#Delta#phi_{cluster-truth} layer_%i", layer), 100, -2e-3, 2e-3);
      h->GetXaxis()->SetTitle("r#Delta#phi_{cluster-truth} (cm)");
      hm->registerHisto(h);
    }

    {
      // rphi cluster errors
      auto h = new TH1F(Form("%srphi_error_%i", get_histo_prefix().c_str(), layer), Form("MVTX r#Delta#phi error layer_%i", layer), 100, 0, 2e-3);
      h->GetXaxis()->SetTitle("r#Delta#phi error (cm)");
      hm->registerHisto(h);
    }

    {
      // phi pulls (cluster - truth)
      auto h = new TH1F(Form("%sphi_pulls_%i", get_histo_prefix().c_str(), layer), Form("MVTX #Delta#phi_{cluster-truth}/#sigma#phi layer_%i", layer), 100, -3, 3);
      h->GetXaxis()->SetTitle("#Delta#phi_{cluster-truth}/#sigma#phi");
      hm->registerHisto(h);
    }

    {
      // z residuals (cluster - truth)
      auto h = new TH1F(Form("%sdz_%i", get_histo_prefix().c_str(), layer), Form("MVTX #Deltaz_{cluster-truth} layer_%i", layer), 100, -3e-3, 3e-3);
      h->GetXaxis()->SetTitle("#Deltaz_{cluster-truth} (cm)");
      hm->registerHisto(h);
    }

    {
      // z cluster errors
      auto h = new TH1F(Form("%sz_error_%i", get_histo_prefix().c_str(), layer), Form("MVTX z error layer_%i", layer), 100, 0, 3e-3);
      h->GetXaxis()->SetTitle("z error (cm)");
      hm->registerHisto(h);
    }

    {
      // z pulls (cluster - truth)
      auto h = new TH1F(Form("%sz_pulls_%i", get_histo_prefix().c_str(), layer), Form("MVTX #Deltaz_{cluster-truth}/#sigmaz layer_%i", layer), 100, -3, 3);
      h->GetXaxis()->SetTitle("#Deltaz_{cluster-truth}/#sigmaz");
      hm->registerHisto(h);
    }

    {
      // total cluster size
      auto h = new TH1F(Form("%sclus_size_%i", get_histo_prefix().c_str(), layer), Form("MVTX cluster size layer_%i", layer), 20, 0, 20);
      h->GetXaxis()->SetTitle("csize");
      hm->registerHisto(h);
    }

    {
      // cluster size in phi
      auto h = new TH1F(Form("%sclus_size_phi_%i", get_histo_prefix().c_str(), layer), Form("MVTX cluster size (#phi) layer_%i", layer), 10, 0, 10);
      h->GetXaxis()->SetTitle("csize_{#phi}");
      hm->registerHisto(h);
    }

    {
      // cluster size in z
      auto h = new TH1F(Form("%sclus_size_z_%i", get_histo_prefix().c_str(), layer), Form("MVTX cluster size (z) layer_%i", layer), 10, 0, 10);
      h->GetXaxis()->SetTitle("csize_{z}");
      hm->registerHisto(h);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int QAG4SimulationMvtx::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK) return res;
  // run evaluation
  evaluate_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
std::string QAG4SimulationMvtx::get_histo_prefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

//________________________________________________________________________
int QAG4SimulationMvtx::load_nodes(PHCompositeNode* topNode)
{
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_cluster_map)
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_cluster_hit_map)
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!m_hit_truth_map)
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_HITTRUTHASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  if (!m_g4hits_mvtx)
  {
    std::cout << PHWHERE << " unable to find DST node G4HIT_MVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
void QAG4SimulationMvtx::evaluate_clusters()
{
  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // load relevant histograms
  struct HistogramList
  {
    TH1* drphi = nullptr;
    TH1* rphi_error = nullptr;
    TH1* phi_pulls = nullptr;

    TH1* dz = nullptr;
    TH1* z_error = nullptr;
    TH1* z_pulls = nullptr;

    TH1* csize = nullptr;
    TH1* csize_phi = nullptr;
    TH1* csize_z = nullptr;
  };

  using HistogramMap = std::map<int, HistogramList>;
  HistogramMap histograms;

  for (const auto& layer : m_layers)
  {
    HistogramList h;
    h.drphi = dynamic_cast<TH1*>(hm->getHisto(Form("%sdrphi_%i", get_histo_prefix().c_str(), layer)));
    h.rphi_error = dynamic_cast<TH1*>(hm->getHisto(Form("%srphi_error_%i", get_histo_prefix().c_str(), layer)));
    h.phi_pulls = dynamic_cast<TH1*>(hm->getHisto(Form("%sphi_pulls_%i", get_histo_prefix().c_str(), layer)));

    h.dz = dynamic_cast<TH1*>(hm->getHisto(Form("%sdz_%i", get_histo_prefix().c_str(), layer)));
    h.z_error = dynamic_cast<TH1*>(hm->getHisto(Form("%sz_error_%i", get_histo_prefix().c_str(), layer)));
    h.z_pulls = dynamic_cast<TH1*>(hm->getHisto(Form("%sz_pulls_%i", get_histo_prefix().c_str(), layer)));

    h.csize = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_%i", get_histo_prefix().c_str(), layer)));
    h.csize_phi = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_phi_%i", get_histo_prefix().c_str(), layer)));
    h.csize_z = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_z_%i", get_histo_prefix().c_str(), layer)));

    histograms.insert(std::make_pair(layer, h));
  }

  for (const auto& hitsetkey : m_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::mvtxId))
  {
    auto range = m_cluster_map->getClusters(hitsetkey);
    for (auto clusterIter = range.first; clusterIter != range.second; ++clusterIter)
    {
      // get cluster key, detector id and check
      const auto& key = clusterIter->first;
      // get cluster
      const auto& cluster = clusterIter->second;
      const auto global = m_tGeometry->getGlobalPosition(key, cluster);

      // get relevant cluster information
      const auto r_cluster = QAG4Util::get_r(global(0), global(1));
      const auto z_cluster = global(2);
      const auto phi_cluster = (float) std::atan2(global(1), global(0));

      double phi_error = 0;
      double z_error = 0;
      if (m_cluster_version == 3)
      {
        phi_error = cluster->getRPhiError() / r_cluster;
        z_error = cluster->getZError();
      }
      else
      {
        auto para_errors = _ClusErrPara.get_si_cluster_error(cluster, key);
        phi_error = sqrt(para_errors.first) / r_cluster;
        z_error = sqrt(para_errors.second);
      }

      // find associated g4hits
      const auto g4hits = find_g4hits(key);

      // get relevant truth information
      const auto x_truth = QAG4Util::interpolate<&PHG4Hit::get_x>(g4hits, r_cluster);
      const auto y_truth = QAG4Util::interpolate<&PHG4Hit::get_y>(g4hits, r_cluster);
      const auto z_truth = QAG4Util::interpolate<&PHG4Hit::get_z>(g4hits, r_cluster);
      const auto phi_truth = std::atan2(y_truth, x_truth);

      const auto dphi = QAG4Util::delta_phi(phi_cluster, phi_truth);
      const auto dz = z_cluster - z_truth;

      // get layer, get histograms
      const auto layer = TrkrDefs::getLayer(key);
      const auto hiter = histograms.find(layer);
      if (hiter == histograms.end()) continue;

      // fill phi residuals, errors and pulls
      auto fill = [](TH1* h, float value)
      { if( h ) h->Fill( value ); };
      fill(hiter->second.drphi, r_cluster * dphi);
      fill(hiter->second.rphi_error, r_cluster * phi_error);
      fill(hiter->second.phi_pulls, dphi / phi_error);

      // fill z residuals, errors and pulls
      fill(hiter->second.dz, dz);
      fill(hiter->second.z_error, z_error);
      fill(hiter->second.z_pulls, dz / z_error);

      // cluster sizes
      // get associated hits
      const auto hit_range = m_cluster_hit_map->getHits(key);
      fill(hiter->second.csize, std::distance(hit_range.first, hit_range.second));

      std::set<int> phibins;
      std::set<int> zbins;
      for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
      {
        const auto& hit_key = hit_iter->second;
        phibins.insert(MvtxDefs::getRow(hit_key));
        zbins.insert(MvtxDefs::getCol(hit_key));
      }

      fill(hiter->second.csize_phi, phibins.size());
      fill(hiter->second.csize_z, zbins.size());
    }
  }
}
//_____________________________________________________________________
QAG4SimulationMvtx::G4HitSet QAG4SimulationMvtx::find_g4hits(TrkrDefs::cluskey cluster_key) const
{
  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

  // loop over hits associated to clusters
  const auto range = m_cluster_hit_map->getHits(cluster_key);
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    // hit key
    const auto& hit_key = iter->second;

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    m_hit_truth_map->getG4Hits(hitset_key, hit_key, g4hit_map);

    // find corresponding g4 hist
    for (auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter)
    {
      // g4hit key
      const auto g4hit_key = truth_iter->second.second;

      // g4 hit
      PHG4Hit* g4hit = (TrkrDefs::getTrkrId(hitset_key) == TrkrDefs::mvtxId) ? m_g4hits_mvtx->findHit(g4hit_key) : nullptr;

      // insert in set
      if (g4hit) out.insert(g4hit);
    }
  }

  return out;
}
