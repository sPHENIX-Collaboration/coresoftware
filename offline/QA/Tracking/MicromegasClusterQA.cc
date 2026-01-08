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
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

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

//_____________________________________________________________________
namespace
{
  //! range adaptor to be able to use range-based for loop
  template <class T>
  class range_adaptor
  {
   public:
    explicit range_adaptor(const T& range)
      : m_range(range)
    {
    }
    const typename T::first_type& begin() { return m_range.first; }
    const typename T::second_type& end() { return m_range.second; }

   private:
    T m_range;
  };

  constexpr int m_max_cluster_count = 10;
  constexpr int m_max_cluster_size = 15;
  constexpr double m_max_cluster_charge = 5e3;

}  // namespace

//____________________________________________________________________________..
MicromegasClusterQA::MicromegasClusterQA(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MicromegasClusterQA::InitRun(PHCompositeNode* topNode)
{
  // print configuration
  std::cout << "MicromegasClusterQA::InitRun - m_use_default_pedestal: " << m_use_default_pedestal << std::endl;
  std::cout << "MicromegasClusterQA::InitRun - m_default_pedestal: " << m_default_pedestal << std::endl;
  std::cout
      << "MicromegasClusterQA::InitRun -"
      << " m_calibration_filename: "
      << (m_calibration_filename.empty() ? "unspecified" : m_calibration_filename)
      << std::endl;

  std::cout << "MicromegasClusterQA::InitRun - m_sample_min: " << m_sample_min << std::endl;
  std::cout << "MicromegasClusterQA::InitRun - m_sample_max: " << m_sample_max << std::endl;

  // read calibrations
  if (!m_calibration_filename.empty())
  {
    m_calibration_data.read(m_calibration_filename);
  }

  // get geometry and keep track of tiles per layer
  auto* geomContainer = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  assert(geomContainer);

  // get number of layers
  m_nlayers = geomContainer->get_NLayers();

  // loop over layers
  bool first = true;
  const auto range = geomContainer->get_begin_end();
  for (const auto& [layer, layergeom] : range_adaptor(range))
  {
    auto* const layergeom_mm = dynamic_cast<CylinderGeomMicromegas*>(layergeom);
    const int ntiles = layergeom_mm->get_tiles_count();
    const auto segmentation = layergeom_mm->get_segmentation_type();

    for (int tile = 0; tile < ntiles; ++tile)
    {
      // generate hitset key get detector name and save
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, tile);
      const auto detector_name = m_mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);
      m_detector_names.push_back(detector_name);
    }

    // keep track of first layer
    if (first)
    {
      first = false;
      m_firstlayer = layer;
    }
  }

  create_histograms();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasClusterQA::process_event(PHCompositeNode* topNode)
{
  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert(m_tGeometry);

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(m_hitsetcontainer);

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert(m_cluster_map);

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  assert(m_cluster_hit_map);

  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // keep track of how many good clusters per detector
  std::array<int, MicromegasDefs::m_nfee> cluster_count = {};
  std::array<int, MicromegasDefs::m_nfee> good_cluster_count = {};

  // first loop over TPOT hitsets
  for (const auto& [hitsetkey, hitset] : range_adaptor(m_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId)))
  {
    // get detector name, layer and tile associated to this hitset key
    const int layer = TrkrDefs::getLayer(hitsetkey);
    const int tile = MicromegasDefs::getTileId(hitsetkey);
    //    const auto detector_name = m_mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

    // detector id
    const int detid = tile + (MicromegasDefs::m_ntiles * (layer - m_firstlayer));

    // get clusters
    const auto cluster_range = m_cluster_map->getClusters(hitsetkey);
    cluster_count[detid] = std::distance(cluster_range.first, cluster_range.second);

    // loop over clusters
    for (const auto& [ckey, cluster] : range_adaptor(cluster_range))
    {
      // find associated hits
      const auto hit_range = m_cluster_hit_map->getHits(ckey);

      // check hit samples
      // if none of the associated hits' sample is within acceptable range, skip the cluster
      if( std::none_of( hit_range.first, hit_range.second,
        [this]( const TrkrClusterHitAssoc::Map::value_type& pair )
        { return MicromegasDefs::getSample( pair.second ) >= m_sample_min &&  MicromegasDefs::getSample( pair.second ) < m_sample_max; } ) )
      { continue; }

      // store cluster size and fill cluster size histogram
      const int cluster_size = std::distance(hit_range.first, hit_range.second);
      m_h_cluster_size->Fill(detid, cluster_size);

      // calculate cluster charge
      double cluster_charge = 0;
      for (const auto& [ckey2, hitkey] : range_adaptor(hit_range))
      {
        // get strip
        const auto strip = MicromegasDefs::getStrip(hitkey);

        // get associated hit
        auto* const hit = hitset->getHit(hitkey);
        assert(hit);

        // get adc, remove pedestal, increment total charge
        const auto adc = hit->getAdc();
        const double pedestal = m_use_default_pedestal ? m_default_pedestal : m_calibration_data.get_pedestal_mapped(hitsetkey, strip);
        cluster_charge += (adc - pedestal);
      }

      // fill cluster charge histogram
      m_h_cluster_charge->Fill(detid, cluster_charge);

      // increment good clusters
      /*
       * we cut on cluster charge > 200 to define good clusters.
       * This is consistent with what has been done offline so far.
       * other cuts considered could include cluster size, cluster strip, etc.
       */
      if (cluster_charge > 200)
      {
        ++good_cluster_count[detid];
      }
    }
  }

  // fill reference and found cluster histograms
  for (int layer = 0; layer < m_nlayers; ++layer)
  {
    for (int tile = 0; tile < MicromegasDefs::m_ntiles; ++tile)
    {
      // get detector id
      const int detid = tile + (MicromegasDefs::m_ntiles * layer);

      // fill multiplicity histogram
      m_h_cluster_multiplicity->Fill(detid, cluster_count[detid]);

      // get reference detector id. It corresponds to the same tile, but on the other layer
      const int detid_ref = detid >= MicromegasDefs::m_ntiles ? detid - MicromegasDefs::m_ntiles : detid + MicromegasDefs::m_ntiles;

      // check if there is exactly one good cluster in the reference detector
      if (good_cluster_count[detid_ref] == 1 && cluster_count[detid_ref] == 1)
      {
        // fill curent detector ref count
        m_h_cluster_count_ref->Fill(detid);

        // if there is one or more cluster in the current detector, also fill good count
        if (cluster_count[detid])
        {
          m_h_cluster_count_found->Fill(detid);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
std::string MicromegasClusterQA::get_histogram_prefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

//____________________________________________________________________________..
void MicromegasClusterQA::create_histograms()
{
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // cluster count histograms
  const int n_detectors = m_detector_names.size();
  m_h_cluster_count_ref = new TH1F(std::format("{}clustercount_ref", get_histogram_prefix()).c_str(), "reference cluster count", n_detectors, 0, n_detectors);
  m_h_cluster_count_found = new TH1F(std::format("{}clustercount_found", get_histogram_prefix()).c_str(), "found cluster count", n_detectors, 0, n_detectors);

  // cluster multiplicity, size and charge distributions
  m_h_cluster_multiplicity = new TH2F(std::format("{}cluster_multiplicity", get_histogram_prefix()).c_str(), "cluster multiplicity", n_detectors, 0, n_detectors, m_max_cluster_count, 0, m_max_cluster_count);
  m_h_cluster_size = new TH2F(std::format("{}cluster_size", get_histogram_prefix()).c_str(), "cluster size", n_detectors, 0, n_detectors, m_max_cluster_size, 0, m_max_cluster_size);
  m_h_cluster_charge = new TH2F(std::format("{}cluster_charge", get_histogram_prefix()).c_str(), "cluster charge", n_detectors, 0, n_detectors, 100, 0, m_max_cluster_charge);

  // assign bin labels
  for (int i = 0; i < n_detectors; ++i)
  {
    for (TH1* h : std::vector<TH1*>({m_h_cluster_count_ref, m_h_cluster_count_found, m_h_cluster_multiplicity, m_h_cluster_size, m_h_cluster_charge}))
    {
      h->GetXaxis()->SetBinLabel(i + 1, m_detector_names[i].c_str());
    }
  }

  // register
  for (TH1* h : std::vector<TH1*>({m_h_cluster_count_ref, m_h_cluster_count_found, m_h_cluster_multiplicity, m_h_cluster_size, m_h_cluster_charge}))
  {
    hm->registerHisto(h);
  }

  return;
}
