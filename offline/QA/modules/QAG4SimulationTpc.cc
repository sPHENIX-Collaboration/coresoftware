#include "QAG4SimulationTpc.h"
#include "QAG4Util.h"
#include "QAHistManagerDef.h"

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4main/PHG4HitContainer.h>

#include <trackbase_historic/ActsTransformations.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>  // for getTrkrId
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4eval/SvtxClusterEval.h>  // for SvtxClusterEval
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTruthEval.h>  // for SvtxTruthEval

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
#include <vector>    // for vector

//________________________________________________________________________
QAG4SimulationTpc::QAG4SimulationTpc(const std::string& name)
  : SubsysReco(name)
  , m_truthContainer(nullptr)
{
}

//________________________________________________________________________
int QAG4SimulationTpc::InitRun(PHCompositeNode* topNode)
{
  // prevent multiple creations of histograms
  if (m_initialized)
    return Fun4AllReturnCodes::EVENT_OK;
  else
    m_initialized = true;

  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(true);
    m_svtxEvalStack->set_verbosity(Verbosity());
  }

  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!m_truthContainer)
  {
    std::cout << "QAG4SimulationTpc::InitRun - Fatal Error - "
              << "unable to find node G4TruthInfo" << std::endl;
    assert(m_truthContainer);
  }

  // find tpc geometry
  PHG4TpcCylinderGeomContainer* geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  // auto geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // TPC has 3 regions, inner, mid and outer
  std::vector<int> region_layer_low = {7, 23, 39};
  ;
  std::vector<int> region_layer_high = {22, 38, 54};

  // get layers from tpc geometry
  // make a layer to region multimap
  const auto range = geom_container->get_begin_end();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    m_layers.insert(iter->first);

    for (int region = 0; region < 3; ++region)
    {
      if (iter->first >= region_layer_low[region] && iter->first <= region_layer_high[region])
        m_layer_region_map.insert(std::make_pair(iter->first, region));
    }
  }

  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create histograms

  // truth clusters
  {
    auto h = new TH1F(Form("%sefficiency_0", get_histo_prefix().c_str()), Form("TPC_truth_clusters"), 48, 7, 54);
    h->GetXaxis()->SetTitle("sPHENIX layer");
    hm->registerHisto(h);
  }
  // matched clusters
  {
    auto h = new TH1F(Form("%sefficiency_1", get_histo_prefix().c_str()), Form("TPC_matched_clusters"), 48, 7, 54);
    h->GetXaxis()->SetTitle("sPHENIX layer");
    hm->registerHisto(h);
  }

  // cluster parameters
  for (int region = 0; region < 3; ++region)
  {
    if (Verbosity()) std::cout << PHWHERE << " adding region " << region << " with layers " << region_layer_low[region] << " to " << region_layer_high[region] << std::endl;
    {
      // rphi residuals (cluster - truth)
      auto h = new TH1F(Form("%sdrphi_%i", get_histo_prefix().c_str(), region), Form("TPC r#Delta#phi_{cluster-truth} region_%i", region), 100, -0.079, 0.075);
      h->GetXaxis()->SetTitle("r#Delta#phi_{cluster-truth} (cm)");
      hm->registerHisto(h);
    }

    {
      // rphi cluster errors
      auto h = new TH1F(Form("%srphi_error_%i", get_histo_prefix().c_str(), region), Form("TPC r#Delta#phi error region_%i", region), 100, 0, 0.075);
      h->GetXaxis()->SetTitle("r#Delta#phi error (cm)");
      hm->registerHisto(h);
    }

    {
      // phi pulls (cluster - truth)
      auto h = new TH1F(Form("%sphi_pulls_%i", get_histo_prefix().c_str(), region), Form("TPC #Delta#phi_{cluster-truth}/#sigma#phi region_%i", region), 100, -3, 3);
      h->GetXaxis()->SetTitle("#Delta#phi_{cluster-truth}/#sigma#phi (cm)");
      hm->registerHisto(h);
    }

    {
      // z residuals (cluster - truth)
      auto h = new TH1F(Form("%sdz_%i", get_histo_prefix().c_str(), region), Form("TPC #Deltaz_{cluster-truth} region_%i", region), 100, -0.19, 0.19);
      h->GetXaxis()->SetTitle("#Delta#z_{cluster-truth} (cm)");
      hm->registerHisto(h);
    }

    {
      // z cluster errors
      auto h = new TH1F(Form("%sz_error_%i", get_histo_prefix().c_str(), region), Form("TPC z error region_%i", region), 100, 0, 0.18);
      h->GetXaxis()->SetTitle("z error (cm)");
      hm->registerHisto(h);
    }

    {
      // z pulls (cluster - truth)
      auto h = new TH1F(Form("%sz_pulls_%i", get_histo_prefix().c_str(), region), Form("TPC #Deltaz_{cluster-truth}/#sigmaz region_%i", region), 100, -3, 3);
      h->GetXaxis()->SetTitle("#Delta#z_{cluster-truth}/#sigmaz (cm)");
      hm->registerHisto(h);
    }

    {
      // total cluster size
      auto h = new TH1F(Form("%sclus_size_%i", get_histo_prefix().c_str(), region), Form("TPC cluster size region_%i", region), 30, 0, 30);
      h->GetXaxis()->SetTitle("csize");
      hm->registerHisto(h);
    }

    {
      // cluster size in phi
      auto h = new TH1F(Form("%sclus_size_phi_%i", get_histo_prefix().c_str(), region), Form("TPC cluster size (#phi) region_%i", region), 10, 0, 10);
      h->GetXaxis()->SetTitle("csize_{#phi}");
      hm->registerHisto(h);
    }

    {
      // cluster size in z
      auto h = new TH1F(Form("%sclus_size_z_%i", get_histo_prefix().c_str(), region), Form("TPC cluster size (z) region_%i", region), 12, 0, 12);
      h->GetXaxis()->SetTitle("csize_{z}");
      hm->registerHisto(h);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int QAG4SimulationTpc::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK) return res;

  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  // run evaluation
  evaluate_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
std::string QAG4SimulationTpc::get_histo_prefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

//________________________________________________________________________
int QAG4SimulationTpc::load_nodes(PHCompositeNode* topNode)
{
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_cluster_map)
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, exiting."
              << std::endl;
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

  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  if (!m_g4hits_tpc)
  {
    std::cout << PHWHERE << " unable to find DST node G4HIT_TPC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
void QAG4SimulationTpc::evaluate_clusters()
{
  SvtxTruthEval* trutheval = m_svtxEvalStack->get_truth_eval();
  assert(trutheval);
  SvtxClusterEval* clustereval = m_svtxEvalStack->get_cluster_eval();
  assert(clustereval);

  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // get histograms for cluster efficiency
  TH1* h_eff0 = dynamic_cast<TH1*>(hm->getHisto(Form("%sefficiency_0", get_histo_prefix().c_str())));
  assert(h_eff0);
  TH1* h_eff1 = dynamic_cast<TH1*>(hm->getHisto(Form("%sefficiency_1", get_histo_prefix().c_str())));
  assert(h_eff1);

  // get histograms for cluster parameters vs truth
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

  for (int region = 0; region < 3; ++region)
  {
    HistogramList h;
    h.drphi = dynamic_cast<TH1*>(hm->getHisto(Form("%sdrphi_%i", get_histo_prefix().c_str(), region)));
    h.rphi_error = dynamic_cast<TH1*>(hm->getHisto(Form("%srphi_error_%i", get_histo_prefix().c_str(), region)));
    h.phi_pulls = dynamic_cast<TH1*>(hm->getHisto(Form("%sphi_pulls_%i", get_histo_prefix().c_str(), region)));

    h.dz = dynamic_cast<TH1*>(hm->getHisto(Form("%sdz_%i", get_histo_prefix().c_str(), region)));
    h.z_error = dynamic_cast<TH1*>(hm->getHisto(Form("%sz_error_%i", get_histo_prefix().c_str(), region)));
    h.z_pulls = dynamic_cast<TH1*>(hm->getHisto(Form("%sz_pulls_%i", get_histo_prefix().c_str(), region)));

    h.csize = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_%i", get_histo_prefix().c_str(), region)));
    h.csize_phi = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_phi_%i", get_histo_prefix().c_str(), region)));
    h.csize_z = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_z_%i", get_histo_prefix().c_str(), region)));

    histograms.insert(std::make_pair(region, h));
  }

  // Get all truth clusters
  //===============
  if (Verbosity() > 0)
    std::cout << PHWHERE << " get all truth clusters for primary particles " << std::endl;

  // PHG4TruthInfoContainer::ConstRange range = m_truthContainer->GetParticleRange();  // all truth cluters
  PHG4TruthInfoContainer::ConstRange range = m_truthContainer->GetPrimaryParticleRange();  // only from primary particles

  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    PHG4Particle* g4particle = iter->second;

    float gtrackID = g4particle->get_track_id();
    float gflavor = g4particle->get_pid();
    float gembed = trutheval->get_embed(g4particle);
    float gprimary = trutheval->is_primary(g4particle);

    if (Verbosity() > 0)
      std::cout << PHWHERE << " PHG4Particle ID " << gtrackID << " gembed " << gembed << " gflavor " << gflavor << " gprimary " << gprimary << std::endl;

    // Get the truth clusters from this particle
    const auto truth_clusters = trutheval->all_truth_clusters(g4particle);

    // get circle fit parameters first
    TrackFitUtils::position_vector_t xy_pts;
    TrackFitUtils::position_vector_t rz_pts;

    for (const auto& [gkey, gclus] : truth_clusters)
    {
      const auto layer = TrkrDefs::getLayer(gkey);
      if (layer < 7) continue;

      float gx = gclus->getX();
      float gy = gclus->getY();
      float gz = gclus->getZ();

      xy_pts.emplace_back(gx, gy);
      rz_pts.emplace_back(std::sqrt(gx * gx + gy * gy), gz);
    }

    // fit a circle through x,y coordinates
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(xy_pts);
    const auto [slope, intercept] = TrackFitUtils::line_fit(rz_pts);

    // skip chain entirely if fit fails
    if (std::isnan(R)) continue;

    // process residuals and pulls
    for (const auto& [gkey, gclus] : truth_clusters)
    {
      const auto layer = TrkrDefs::getLayer(gkey);
      const auto detID = TrkrDefs::getTrkrId(gkey);
      // if (detID != TrkrDefs::tpcId) continue;
      if (layer < 7) continue;

      float gx = gclus->getX();
      float gy = gclus->getY();
      float gz = gclus->getZ();
      float gedep = gclus->getError(0, 0);
      float ng4hits = gclus->getAdc();

      const auto gr = QAG4Util::get_r(gclus->getX(), gclus->getY());
      const auto gphi = std::atan2(gclus->getY(), gclus->getX());

      if (Verbosity() > 0)
      {
        std::cout << "     gkey " << gkey << " detID " << detID << " tpcId " << TrkrDefs::tpcId << " layer " << layer << "  truth clus " << gkey << " ng4hits " << ng4hits << " gr " << gr << " gx " << gx << " gy " << gy << " gz " << gz
                  << " gphi " << gphi << " gedep " << gedep << std::endl;
      }

      // fill the truth cluster histo
      h_eff0->Fill(layer);

      // find matching reco cluster histo
      const auto [rkey, rclus] = clustereval->reco_cluster_from_truth_cluster(gkey, gclus);
      if (rclus)
      {
        // fill the matched cluster histo
        h_eff1->Fill(layer);

        const auto global = m_tGeometry->getGlobalPosition(rkey, rclus);

        // get relevant cluster information
        const auto r_cluster = QAG4Util::get_r(global(0), global(1));
        const auto z_cluster = global(2);
        const auto phi_cluster = (float) std::atan2(global(1), global(0));

        double phi_error = 0;
        double z_error = 0;

        if (m_cluster_version == 3)
        {
          phi_error = rclus->getRPhiError() / r_cluster;
          z_error = rclus->getZError();
        }
        else
        {
          float r = r_cluster;
          double alpha = (r * r) / (2 * r * R);
          double beta = slope;

          auto para_errors = _ClusErrPara.get_cluster_error(rclus, rkey, alpha, beta);
          phi_error = sqrt(para_errors.first) / r_cluster;
          z_error = sqrt(para_errors.second);
        }

        const auto dphi = QAG4Util::delta_phi(phi_cluster, gphi);
        const auto dz = z_cluster - gz;

        // get region from layer, fill histograms
        const auto it = m_layer_region_map.find(layer);
        int region = it->second;

        if (Verbosity() > 0)
        {
          std::cout << "   Found match in layer " << layer << " region " << region << " for gtrackID " << gtrackID << std::endl;
          std::cout << "      x " << rclus->getX() << " y " << rclus->getY() << " z " << rclus->getZ() << std::endl;
          std::cout << "     gx " << gclus->getX() << " gy " << gclus->getY() << " gz " << gclus->getZ() << std::endl;
          std::cout << "     drphi " << r_cluster * dphi << " rphi_error " << r_cluster * phi_error << " dz " << dz << " z_error " << z_error << std::endl;
        }

        const auto hiter = histograms.find(region);
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
        const auto hit_range = m_cluster_hit_map->getHits(rkey);
        fill(hiter->second.csize, std::distance(hit_range.first, hit_range.second));

        std::set<int> phibins;
        std::set<int> zbins;
        for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
        {
          const auto& hit_key = hit_iter->second;
          phibins.insert(TpcDefs::getPad(hit_key));
          zbins.insert(TpcDefs::getTBin(hit_key));
        }

        fill(hiter->second.csize_phi, phibins.size());
        fill(hiter->second.csize_z, zbins.size());
      }
    }
  }
}

//_____________________________________________________________________
QAG4SimulationTpc::G4HitSet QAG4SimulationTpc::find_g4hits(TrkrDefs::cluskey cluster_key) const
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
      PHG4Hit* g4hit = (TrkrDefs::getTrkrId(hitset_key) == TrkrDefs::tpcId) ? m_g4hits_tpc->findHit(g4hit_key) : nullptr;

      // insert in set
      if (g4hit) out.insert(g4hit);
    }
  }

  return out;
}
