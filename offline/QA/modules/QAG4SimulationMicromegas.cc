#include "QAG4SimulationMicromegas.h"
#include "QAG4Util.h"
#include "QAHistManagerDef.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <mvtx/MvtxDefs.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>  // for getTrkrId, getHit...
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <micromegas/CylinderGeomMicromegas.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TAxis.h>  // for TAxis
#include <TH1.h>
#include <TString.h>  // for Form
#include <TVector3.h>

#include <cassert>
#include <iostream>  // for operator<<, basic...
#include <iterator>  // for distance
#include <map>       // for map
#include <utility>   // for pair, make_pair

//_____________________________________________________________________
namespace 
{
 
  //! square
  template<class T> inline constexpr T square( T x ) { return x*x; }
  
  //! needed for weighted linear interpolation
  struct interpolation_data_t 
  {
    using list = std::vector<interpolation_data_t>;
    double x() const { return position.x(); }
    double y() const { return position.y(); }
    double z() const { return position.z(); }
    
    double px() const { return momentum.x(); }
    double py() const { return momentum.y(); }
    double pz() const { return momentum.z(); }
    
    TVector3 position;
    TVector3 momentum;
    double weight = 1;
  };
  
 //! rotate size or covariant matrix to polar coordinates and return component along provided phi angle
 /** this is copied from TrkrClusterv2.cc, with the difference that one must use the tile azimuthal angle instead of that of the cluster */
 template<float (TrkrCluster::*accessor)(unsigned int, unsigned int) const>
    float rotate( const double phi, const TrkrCluster* cluster )
  {
    const auto cosphi = std::cos(phi);
    const auto sinphi = std::sin(phi);

    return
      square(sinphi)*(cluster->*accessor)(0,0) +
      square(cosphi)*(cluster->*accessor)(1,1) +
      2.*cosphi*sinphi*(cluster->*accessor)(0,1);
  }
  
  //! calculate the interpolation of member function called on all members in collection to the provided y_extrap
  template<double (interpolation_data_t::*accessor)() const>
  double interpolate( const interpolation_data_t::list& hits, double y_extrap )
  {
    
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swy = 0;
    double swy2 = 0;
    double swx = 0;
    double swyx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x = (hit.*accessor)();
      const double w = hit.weight;
      if( w <= 0 ) continue;

      valid = true;
      const double y = hit.y();
        
      sw += w;
      swy += w*y;
      swy2 += w*square(y);
      swx += w*x;
      swyx += w*x*y;
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swyx - swy*swx);
    const auto beta = (swy2*swx - swy*swyx);
    const auto denom = (sw*swy2 - square(swy));

    return ( alpha*y_extrap + beta )/denom;
  }
  
}

//________________________________________________________________________
QAG4SimulationMicromegas::QAG4SimulationMicromegas(const std::string& name)
  : SubsysReco(name)
{
}

//________________________________________________________________________
int QAG4SimulationMicromegas::InitRun(PHCompositeNode* topNode)
{
  // prevent multiple creations of histograms
  if (m_initialized)
    return Fun4AllReturnCodes::EVENT_OK;
  else
    m_initialized = true;

  // find mvtx geometry
  auto geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!geom_container)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERGEOM_MICROMEGAS_FULL" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get layers from mvtx geometry
  const auto range = geom_container->get_begin_end();
  for (auto iter = range.first; iter != range.second; ++iter)
  { m_layers.insert(iter->first); }

  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create histograms
  for (const auto& layer : m_layers)
  {
    if (Verbosity()) std::cout << PHWHERE << " adding layer " << layer << std::endl;

    // get layer geometry
    const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(geom_container->GetLayerGeom(layer));
    assert( layergeom );

    // get segmentation type
    const auto segmentation_type = layergeom->get_segmentation_type();
    const bool is_segmentation_phi = (segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI);
    {
      // residuals (cluster - truth)
      const double max_residual = is_segmentation_phi ? 0.04:0.08;
      auto h = new TH1F(Form("%sresidual_%i", get_histo_prefix().c_str(), layer),
        Form("micromegas %s_{cluster-truth} layer_%i",
        is_segmentation_phi ? "r#Delta#phi":"#Deltaz", layer), 100, -max_residual, max_residual );
      h->GetXaxis()->SetTitle( Form( "%s_{cluster-truth} (cm)", is_segmentation_phi ? "r#Delta#phi":"#Deltaz" ) );
      hm->registerHisto(h);
    }

    {
      // cluster errors
      const double max_error =  is_segmentation_phi ? 0.04:0.08;
      auto h = new TH1F(Form("%sresidual_error_%i", get_histo_prefix().c_str(), layer),
        Form("micromegas %s error layer_%i",
        is_segmentation_phi ? "r#Delta#phi":"#Deltaz", layer), 100, 0, max_error);
      h->GetXaxis()->SetTitle( Form( "%s error (cm)", is_segmentation_phi ? "r#Delta#phi":"#Deltaz" ) );
      hm->registerHisto(h);
    }

    {
      // pulls (cluster - truth)
      auto h = new TH1F(Form("%scluster_pulls_%i", get_histo_prefix().c_str(), layer),
        Form("micromegas %s layer_%i",
        is_segmentation_phi ? "#Delta#phi/#sigma#phi":"#Deltaz/#sigmaz", layer), 100, -5, 5);
      h->GetXaxis()->SetTitle(is_segmentation_phi ? "#Delta#phi/#sigma#phi":"#Deltaz/#sigmaz");
      hm->registerHisto(h);
    }

    {
      // cluster size
      auto h = new TH1F(Form("%sclus_size_%i", get_histo_prefix().c_str(), layer), Form("micromegas cluster size layer_%i", layer), 20, 0, 20);
      h->GetXaxis()->SetTitle("csize");
      hm->registerHisto(h);
    }

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int QAG4SimulationMicromegas::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK) return res;
  // run evaluation
  evaluate_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
std::string QAG4SimulationMicromegas::get_histo_prefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

//________________________________________________________________________
int QAG4SimulationMicromegas::load_nodes(PHCompositeNode* topNode)
{

  m_micromegas_geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL" );
  if (!m_micromegas_geonode)
  {
    std::cout << PHWHERE << " ERROR: Can't find CYLINDERGEOM_MICROMEGAS_FULL." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hitsets)
  {
    std::cout << PHWHERE << " ERROR: Can't find TrkrHitSetContainer." << std::endl;
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

  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");
  if (!m_g4hits_micromegas)
  {
    std::cout << PHWHERE << " unable to find DST node G4HIT_MVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
void QAG4SimulationMicromegas::evaluate_clusters()
{
  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // load relevant histograms
  struct HistogramList
  {
    TH1* residual = nullptr;
    TH1* residual_error = nullptr;
    TH1* pulls = nullptr;

    TH1* csize = nullptr;
  };

  using HistogramMap = std::map<int, HistogramList>;
  HistogramMap histograms;

  for (const auto& layer : m_layers)
  {
    HistogramList h;
    h.residual = dynamic_cast<TH1*>(hm->getHisto(Form("%sresidual_%i", get_histo_prefix().c_str(), layer)));
    h.residual_error = dynamic_cast<TH1*>(hm->getHisto(Form("%sresidual_error_%i", get_histo_prefix().c_str(), layer)));
    h.pulls = dynamic_cast<TH1*>(hm->getHisto(Form("%scluster_pulls_%i", get_histo_prefix().c_str(), layer)));
    h.csize = dynamic_cast<TH1*>(hm->getHisto(Form("%sclus_size_%i", get_histo_prefix().c_str(), layer)));

    histograms.insert(std::make_pair(layer, h));
  }

  // loop over hitsets
  const auto hitsetrange = m_hitsets->getHitSets(TrkrDefs::TrkrId::micromegasId);
  for (auto hitsetiter = hitsetrange.first; hitsetiter != hitsetrange.second; ++hitsetiter)
  {

    // get hitsetkey, layer and tileid
    const auto hitsetkey = hitsetiter->first;
    
    // get layer
    const auto layer = TrkrDefs::getLayer(hitsetkey);

    // get tileid
    const auto tileid = MicromegasDefs::getTileId(hitsetkey);
    
    // load geometry
    const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geonode->GetLayerGeom(layer));
    if( !layergeom ) continue;
    
    // get relevant histograms
    const auto hiter = histograms.find(layer);
    if (hiter == histograms.end()) continue;
    
    // get associated clusters
    const auto range = m_cluster_map->getClusters(hitsetkey);
    for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter ){

      // get cluster key
      const auto& key = clusterIter->first;

      // get cluster
      const auto& cluster = clusterIter->second;

      // get segmentation type
      const auto segmentation_type = MicromegasDefs::getSegmentationType( key );

      // get relevant cluster information
      const auto phi = -layergeom->get_center_phi( tileid );
      const auto rphi_error = std::sqrt(rotate<&TrkrCluster::getError>(phi, cluster) );        
      const auto z_error = cluster->getZError();

      // convert cluster position to local tile coordinates
      const TVector3 cluster_world( cluster->getX(), cluster->getY(), cluster->getZ() );
      const TVector3 cluster_local = layergeom->get_local_from_world_coords( tileid, cluster_world );

      // find associated g4hits
      const auto g4hits = find_g4hits(key);
 
      // convert hits to list of interpolation_data_t
      interpolation_data_t::list hits;
      for( const auto& g4hit:g4hits )
      {
        const auto weight = g4hit->get_edep();
        for( int i = 0; i < 2; ++i )
        {
          
          // convert position to local
          TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
          TVector3 g4hit_local = layergeom->get_local_from_world_coords( tileid, g4hit_world );
          
          // convert momentum to local
          TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
          TVector3 momentum_local = layergeom->get_local_from_world_vect( tileid, momentum_world );
          
          hits.push_back( {.position = g4hit_local, .momentum = momentum_local, .weight = weight } );
        }
      }

      // do position interpolation
      const auto y_extrap = cluster_local.y();
      const TVector3 interpolation_local( 
        interpolate<&interpolation_data_t::x>( hits, y_extrap ),
        interpolate<&interpolation_data_t::y>( hits, y_extrap ),
        interpolate<&interpolation_data_t::z>( hits, y_extrap ) );

      // fill phi residuals, errors and pulls
      auto fill = [](TH1* h, float value) { if( h ) h->Fill( value ); };
      switch( segmentation_type )
      {
        case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
        {
          const auto drphi = cluster_local.x() - interpolation_local.x();
          fill(hiter->second.residual, drphi );
          fill(hiter->second.residual_error, rphi_error);
          fill(hiter->second.pulls, drphi / rphi_error );
          break;
        }
        
        case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
        {
          const auto dz = cluster_local.z() - interpolation_local.z();
          fill(hiter->second.residual, dz);
          fill(hiter->second.residual_error, z_error);
          fill(hiter->second.pulls, dz / z_error);
          break;
        }
      }
      
      // cluster size
      // get associated hits
      const auto hit_range = m_cluster_hit_map->getHits(key);
      fill(hiter->second.csize, std::distance(hit_range.first, hit_range.second));

    }
  }
}
//_____________________________________________________________________
QAG4SimulationMicromegas::G4HitSet QAG4SimulationMicromegas::find_g4hits(TrkrDefs::cluskey cluster_key) const
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
      PHG4Hit* g4hit = (TrkrDefs::getTrkrId(hitset_key) == TrkrDefs::micromegasId) ? m_g4hits_micromegas->findHit(g4hit_key) : nullptr;

      // insert in set
      if (g4hit) out.insert(g4hit);
    }
  }

  return out;
}
