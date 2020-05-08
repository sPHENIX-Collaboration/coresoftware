#include "QAG4SimulationIntt.h"
#include "QAG4Util.h"
#include "QAHistManagerDef.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <TH1F.h>

#include <cassert>

//________________________________________________________________________
QAG4SimulationIntt::QAG4SimulationIntt(const std::string &name)
  : SubsysReco(name)
{}

//________________________________________________________________________
int QAG4SimulationIntt::InitRun(PHCompositeNode *topNode)
{

  // prevent multiple creations of histograms
  if( m_initialized ) return Fun4AllReturnCodes::EVENT_OK;
  else m_initialized = true;

  // find intt geometry
  auto geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container)
  {
    std::cout << PHWHERE << " unable to find DST node CYLINDERGEOM_INTT" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get layers from intt geometry
  const auto range = geom_container->get_begin_end();
  for( auto iter = range.first; iter != range.second; ++iter )
  { m_layers.insert( iter->first ); }
  
  // histogram manager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create histograms
  for( const auto& layer:m_layers )
  {
    std::cout << PHWHERE << " adding layer " << layer << std::endl;
    {
      // rphi residuals (cluster - truth)
      auto h = new TH1F( Form( "%sdrphi_%i", get_histo_prefix().c_str(), layer ), Form( "INTT r#Delta#phi_{cluster-truth} layer_%i", layer ), 100, -1e-2, 1e-2 );
      h->GetXaxis()->SetTitle( "r#Delta#phi_{cluster-truth} (cm)" );
      hm->registerHisto(h);
    }

    {
      // rphi cluster errors
      auto h = new TH1F( Form( "%srphi_error_%i", get_histo_prefix().c_str(), layer ), Form( "INTT r#Delta#phi error layer_%i", layer ), 100, 0, 1e-2 );
      h->GetXaxis()->SetTitle( "r#Delta#phi error (cm)" );
      hm->registerHisto(h);
    }

    {
      // phi pulls (cluster - truth)
      auto h = new TH1F( Form( "%sphi_pulls_%i", get_histo_prefix().c_str(), layer ), Form( "INTT #Delta#phi_{cluster-truth}/#sigma#phi layer_%i", layer ), 100, -3, 3 );
      h->GetXaxis()->SetTitle( "#Delta#phi_{cluster-truth}/#sigma#phi (cm)" );
      hm->registerHisto(h);
    }

    {
      // z residuals (cluster - truth)
      auto h = new TH1F( Form( "%sdz_%i", get_histo_prefix().c_str(), layer ), Form( "INTT #Deltaz_{cluster-truth} layer_%i", layer ), 100, -2.5, 2.5 );
      h->GetXaxis()->SetTitle( "#Delta#z_{cluster-truth} (cm)" );
      hm->registerHisto(h);
    }

    {
      // z cluster errors
      auto h = new TH1F( Form( "%sz_error_%i", get_histo_prefix().c_str(), layer ), Form( "INTT z error layer_%i", layer ), 100, 0, 2.5 );
      h->GetXaxis()->SetTitle( "z error (cm)" );
      hm->registerHisto(h);
    }

    {
      // z pulls (cluster - truth)
      auto h = new TH1F( Form( "%sz_pulls_%i", get_histo_prefix().c_str(), layer ), Form( "INTT #Deltaz_{cluster-truth}/#sigmaz layer_%i", layer ), 100, -3, 3 );
      h->GetXaxis()->SetTitle( "#Delta#z_{cluster-truth}/#sigmaz (cm)" );
      hm->registerHisto(h);
    }

  }

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int QAG4SimulationIntt::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // run evaluation
  evaluate_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
std::string QAG4SimulationIntt::get_histo_prefix() const
{ return std::string("h_") + Name() + std::string("_"); }

//________________________________________________________________________
int QAG4SimulationIntt::load_nodes( PHCompositeNode* topNode )
{
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if( !m_cluster_map )
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if( !m_cluster_hit_map )
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  if( !m_hit_truth_map )
  {
    std::cout << PHWHERE << " unable to find DST node TRKR_HITTRUTHASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  if( !m_g4hits_intt )
  {
    std::cout << PHWHERE << " unable to find DST node G4HIT_INTT" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//________________________________________________________________________
void QAG4SimulationIntt::evaluate_clusters()
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
  };

  using HistogramMap = std::map<int,HistogramList>;
  HistogramMap histograms;

  for( const auto& layer: m_layers )
  {
    HistogramList h;
    h.drphi = dynamic_cast<TH1*>( hm->getHisto(Form( "%sdrphi_%i", get_histo_prefix().c_str(), layer )) );
    h.rphi_error = dynamic_cast<TH1*>( hm->getHisto(Form( "%srphi_error_%i", get_histo_prefix().c_str(), layer )) );
    h.phi_pulls = dynamic_cast<TH1*>( hm->getHisto(Form( "%sphi_pulls_%i", get_histo_prefix().c_str(), layer )) );

    h.dz = dynamic_cast<TH1*>( hm->getHisto(Form( "%sdz_%i", get_histo_prefix().c_str(), layer )) );
    h.z_error = dynamic_cast<TH1*>( hm->getHisto(Form( "%sz_error_%i", get_histo_prefix().c_str(), layer )) );
    h.z_pulls = dynamic_cast<TH1*>( hm->getHisto(Form( "%sz_pulls_%i", get_histo_prefix().c_str(), layer )) );

    histograms.insert( std::make_pair(layer,h));
  }

  auto range = m_cluster_map->getClusters();
  for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
  {

    // get cluster key, detector id and check
    const auto& key = clusterIter->first;
    const auto detId = TrkrDefs::getTrkrId( key );
    if( detId != TrkrDefs::inttId ) continue;

    // get cluster
    const auto& cluster = clusterIter->second;

    // get relevant cluster information
    const auto r_cluster = QAG4Util::get_r( cluster->getX(), cluster->getY() );
    const auto z_cluster = cluster->getZ();
    const auto phi_cluster = std::atan2( cluster->getY(), cluster->getX() );
    const auto phi_error = cluster->getPhiError();
    const auto z_error = cluster->getZError();

    // find associated g4hits
    const auto g4hits = find_g4hits( key );

    // get relevant truth information
    const auto x_truth = QAG4Util::interpolate<&PHG4Hit::get_x>( g4hits, r_cluster );
    const auto y_truth = QAG4Util::interpolate<&PHG4Hit::get_y>( g4hits, r_cluster );
    const auto z_truth = QAG4Util::interpolate<&PHG4Hit::get_z>( g4hits, r_cluster );
    const auto phi_truth = std::atan2( y_truth, x_truth );

    const auto dphi = QAG4Util::delta_phi( phi_cluster, phi_truth );
    const auto dz = z_cluster-z_truth;

    // get layer, get histograms
    const auto layer = TrkrDefs::getLayer(key);
    const auto hiter = histograms.find( layer );
    if( hiter == histograms.end() ) continue;

    // fill histograms
    auto fill = []( TH1* h, float value ) { if( h ) h->Fill( value ); };
    fill( hiter->second.drphi, r_cluster*dphi );
    fill( hiter->second.rphi_error, r_cluster*phi_error );
    fill( hiter->second.phi_pulls, dphi/phi_error );

    fill( hiter->second.dz, dz );
    fill( hiter->second.z_error, z_error );
    fill( hiter->second.z_pulls, dz/z_error );

  }
}

//_____________________________________________________________________
QAG4SimulationIntt::G4HitSet QAG4SimulationIntt::find_g4hits( TrkrDefs::cluskey cluster_key ) const
{

  // check maps
  if( !( m_cluster_hit_map && m_hit_truth_map ) ) return G4HitSet();

  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

  // loop over hits associated to clusters
  const auto range = m_cluster_hit_map->getHits(cluster_key);
  for( auto iter = range.first; iter != range.second; ++iter )
  {

    // hit key
    const auto& hit_key = iter->second;

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    m_hit_truth_map->getG4Hits( hitset_key, hit_key, g4hit_map );

    // find corresponding g4 hist
    for( auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter )
    {

      // g4hit key
      const auto g4hit_key = truth_iter->second.second;

      // g4 hit
      PHG4Hit* g4hit = ( TrkrDefs::getTrkrId( hitset_key ) == TrkrDefs::inttId ) ?
        m_g4hits_intt->findHit( g4hit_key ):nullptr;

      // insert in set
      if( g4hit ) out.insert( g4hit );

    }
  }

  return out;
}
