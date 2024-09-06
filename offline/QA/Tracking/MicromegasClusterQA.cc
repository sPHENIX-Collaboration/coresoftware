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

#include <boost/format.hpp>


//_____________________________________________________________________
namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    explicit range_adaptor( const T& range ):m_range(range){}
    const typename T::first_type& begin() {return m_range.first;}
    const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };
}

//____________________________________________________________________________..
MicromegasClusterQA::MicromegasClusterQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MicromegasClusterQA::InitRun(PHCompositeNode *topNode)
{
  // print configuration
  std::cout << "MicromegasClusterQA::InitRun - m_use_default_pedestal: " << m_use_default_pedestal << std::endl;
  std::cout << "MicromegasClusterQA::InitRun - m_default_pedestal: " << m_default_pedestal << std::endl;
  std::cout
    << "MicromegasClusterQA::InitRun -"
    << " m_calibration_filename: "
    << (m_calibration_filename.empty() ? "unspecified":m_calibration_filename )
    << std::endl;

  // read calibrations
  if( !m_calibration_filename.empty() )
  { m_calibration_data.read( m_calibration_filename ); }

  // get geometry and keep track of tiles per layer
  auto geomContainer = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  assert( geomContainer );

  // get number of layers
  m_nlayers = geomContainer->get_NLayers();

  // loop over layers
  bool first = true;
  const auto range = geomContainer->get_begin_end();
  for( const auto& [layer, layergeom]:range_adaptor( range ) )
  {
    const auto layergeom_mm = static_cast<CylinderGeomMicromegas*>(layergeom);
    const int ntiles = layergeom_mm->get_tiles_count();
    const auto segmentation = layergeom_mm->get_segmentation_type();

    for( int tile=0; tile<ntiles; ++tile )
    {
      // generate hitset key get detector name and save
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, tile );
      const auto detector_name = m_mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      m_detector_names.push_back(detector_name);
    }

    // keep track of first layer
    if( first )
    {
      first=false;
      m_firstlayer = layer;
    }
  }

  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasClusterQA::process_event(PHCompositeNode *topNode)
{

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(m_hitsetcontainer);

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert( m_cluster_map );

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  assert( m_cluster_hit_map );

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // keep track of how many good clusters per detector
  std::array<int, MicromegasDefs::m_nfee> cluster_count = {};
  std::array<int, MicromegasDefs::m_nfee> good_cluster_count = {};

  // first loop over TPOT hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId)))
  {
    // get detector name, layer and tile associated to this hitset key
    const int layer = TrkrDefs::getLayer(hitsetkey);
    const int tile = MicromegasDefs::getTileId(hitsetkey);
    const auto detector_name = m_mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);

    // detector id
    const int detid = tile + MicromegasDefs::m_ntiles*(layer-m_firstlayer);

    // get profile histogram
    const auto h_profile = dynamic_cast<TH2 *>(hm->getHisto((boost::format("%sncluspertile%i_%i") % getHistoPrefix()%(layer-m_firstlayer)%tile).str()));

    // get clusters
    const auto cluster_range = m_cluster_map->getClusters(hitsetkey);
    cluster_count[detid] = std::distance( cluster_range.first, cluster_range.second );

    // loop over clusters
    for( const auto& [ckey,cluster]:range_adaptor(cluster_range))
    {
      // fill profile
      h_profile->Fill(cluster->getLocalY(), cluster->getLocalX());

      // find associated hits
      const auto hit_range = m_cluster_hit_map->getHits(ckey);

      // store cluster size
      // const int cluster_size = std::distance( hit_range.first, hit_range.second );

      // calculate cluster charge
      double cluster_charge = 0;
      for( const auto& [ckey2, hitkey]:range_adaptor( hit_range ) )
      {

        // get strip
        const auto strip = MicromegasDefs::getStrip( hitkey );

        // get associated hit
        const auto hit = hitset->getHit( hitkey );
        assert( hit );

        // get adc, remove pedestal, increment total charge
        const auto adc = hit->getAdc();
        const double pedestal = m_use_default_pedestal ?
          m_default_pedestal:
          m_calibration_data.get_pedestal_mapped(hitsetkey, strip);
        cluster_charge += (adc-pedestal);
      }

      // increment good clusters
      /*
       * we cut on cluster charge > 200 to define good clusters.
       * This is consistent with what has been done offline so far.
       * other cuts considered could include cluster size, cluster strip, etc.
       */
      if( cluster_charge > 200 )
      { ++good_cluster_count[detid]; }
    }
  }

  // fill reference and found cluster histograms
  for( int layer = 0; layer < m_nlayers; ++layer )
  {
    for( int tile =0; tile < MicromegasDefs::m_ntiles; ++tile )
    {
      // get detector id
      const int detid = tile+MicromegasDefs::m_ntiles*layer;

      // get reference detector id. It corresponds to the same tile, but on the other layer
      const int detid_ref = detid >= MicromegasDefs::m_ntiles ?
        detid-MicromegasDefs::m_ntiles:
        detid+MicromegasDefs::m_ntiles;

      // check if there is exactly one good cluster in the reference detector
      if(good_cluster_count[detid_ref]==1 && cluster_count[detid_ref]==1)
      {
        // fill curent detector ref count
        m_h_clustercount_ref->Fill(detid);

        // if there is one or more cluster in the current detector, also fill good count
        if(cluster_count[detid])
        { m_h_clustercount_found->Fill(detid); }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasClusterQA::EndRun(const int /*runnumber*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//____________________________________________________________________________..
std::string MicromegasClusterQA::getHistoPrefix() const
{ return std::string("h_") + Name() + std::string("_"); }

//____________________________________________________________________________..
void MicromegasClusterQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // cluster count histograms
  const int n_detectors = m_detector_names.size();
  m_h_clustercount_ref = new TH1F( (boost::format("%sclustercount_ref") % getHistoPrefix()).str().c_str(), "reference cluster count", n_detectors, 0, n_detectors );
  m_h_clustercount_found = new TH1F( (boost::format("%sclustercount_found") % getHistoPrefix()).str().c_str(), "found cluster count", n_detectors, 0, n_detectors );

  // cluster positions vs tile
  for( int i = 0; i < n_detectors; ++i )
  {
    m_h_clustercount_ref->GetXaxis()->SetBinLabel(i+1, m_detector_names[i].c_str());
    m_h_clustercount_found->GetXaxis()->SetBinLabel(i+1, m_detector_names[i].c_str());
  }

  // register
  hm->registerHisto(m_h_clustercount_ref);
  hm->registerHisto(m_h_clustercount_found);

  // cluster positions vs tile
  for( int layer = 0; layer < m_nlayers; ++layer )
  {
    for (int tile = 0; tile < MicromegasDefs::m_ntiles; tile++)
    {
      auto h = new TH2F(
        (boost::format("%sncluspertile%i_%i") % getHistoPrefix() % layer% tile).str().c_str(),
        "Micromegas clusters per tile", 2000, -30, 30, 2000, -20, 20);
      h->GetXaxis()->SetTitle("Local x [cm]");
      h->GetYaxis()->SetTitle("Local y [cm]");
      hm->registerHisto(h);
    }
  }

  return;
}
