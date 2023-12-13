#include "MicromegasClusterEvaluator_hp.h"

#include "MicromegasGeometryContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasMapping.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/PHRandomSeed.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <gsl/gsl_randist.h>

//_____________________________________________________________________
namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    const typename T::first_type& begin() {return m_range.first;}
    const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };

  //! square
  template<class T> inline constexpr T square( T x ) { return x*x; }

  //* converninece trait for underlying type
  template<class T>
    using underlying_type_t = typename std::underlying_type<T>::type;

  //* convert an strong type enum to integral type
  template<class T>
    constexpr underlying_type_t<T>
    to_underlying_type(T value) noexcept
  { return static_cast<underlying_type_t<T>>(value);}

  // TVector3 streamer
  inline std::ostream& operator << (std::ostream& out, const TVector3& v )
  {
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }

}

//_____________________________________________________________________
void MicromegasClusterEvaluator_hp::Container::Reset()
{
  n_waveforms_signal = 0;
  n_clusters = 0;
  n_phi_clusters = 0;
  n_z_clusters = 0;
  clusters.clear();
  n_detector_clusters.clear();
  n_region_clusters.clear();
  min_cluster_size.clear();
  min_cluster_charge.clear();
  first_cluster_strip.clear();
}

//_____________________________________________________________________
MicromegasClusterEvaluator_hp::MicromegasClusterEvaluator_hp( const std::string& name ):
  SubsysReco( name)
  {}

//_____________________________________________________________________
int MicromegasClusterEvaluator_hp::Init(PHCompositeNode* topNode )
{
  // print configuration
  std::cout << "MicromegasClusterEvaluator_hp::Init - m_use_default_pedestal: " << m_use_default_pedestal << std::endl;
  std::cout << "MicromegasClusterEvaluator_hp::Init - m_default_pedestal: " << m_default_pedestal << std::endl;
  std::cout
    << "MicromegasClusterEvaluator_hp::Init -"
    << " m_calibration_filename: "
    << (m_calibration_filename.empty() ? "unspecified":m_calibration_filename )
    << std::endl;

  // read calibrations
  if( !m_calibration_filename.empty() )
  { m_calibration_data.read( m_calibration_filename ); }

  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MicromegasClusterEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "MicromegasClusterEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  auto newNode = new PHIODataNode<PHObject>( new Container, "MicromegasClusterEvaluator_hp::Container","PHObject");
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasClusterEvaluator_hp::InitRun(PHCompositeNode* topNode )
{
  // load geometry
  PHG4CylinderGeomContainer* geonode = nullptr;
  for( const auto& geonodename: {"CYLINDERGEOM_MICROMEGAS_FULL", "CYLINDERGEOM_MICROMEGAS" } )
  { if(( geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename))) break; }
  assert(geonode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasClusterEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;
  if( m_container ) m_container->Reset();

  evaluate_clusters();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasClusterEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int MicromegasClusterEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // geometry
  const std::string geonodename = "CYLINDERGEOM_MICROMEGAS_FULL";
  m_geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  assert(m_geonode);

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(m_hitsetcontainer);

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if( !m_cluster_map )
  { m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER"); }
  assert( m_cluster_map );

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  assert( m_cluster_hit_map );

  // local container
  m_container = findNode::getClass<Container>(topNode, "MicromegasClusterEvaluator_hp::Container");
  assert(m_container);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void MicromegasClusterEvaluator_hp::evaluate_clusters()
{
  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;

  // clear array
  m_container->Reset();

  m_container->n_detector_clusters.resize(16,0);
  m_container->n_region_clusters.resize(64,0);
  m_container->min_cluster_size.resize(16,0);
  m_container->min_cluster_charge.resize(16,0);
  m_container->first_cluster_strip.resize(16,0);

  // first loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    // process only micromegas clusters
    const bool is_micromegas( TrkrDefs::getTrkrId(hitsetkey) == TrkrDefs::micromegasId );
    if( !is_micromegas ) continue;

    // increment total number of hits
    m_container->n_waveforms_signal += hitset->size();

    // smallest cluster charge in detector
    unsigned short min_cluster_size = 0;

    // smallest cluster charge in detector
    double min_cluster_charge = 0;

    const auto layer = TrkrDefs::getLayer(hitsetkey);
    const auto tile = MicromegasDefs::getTileId(hitsetkey);
    const int detid = int(tile) + 8*(int(layer)-55);

    const auto cluster_range = m_cluster_map->getClusters(hitsetkey);
    const auto nclusters = std::distance( cluster_range.first, cluster_range.second );

    m_container->n_clusters += nclusters;

    if( Verbosity() )
    {
      std::cout
        << "MicromegasClusterEvaluator_hp::evaluate_clusters -"
        << " layer: " << int(layer) << " tile: " << int(tile) << " detid: " << detid
        << " hits: " << hitset->size()
        << " clusters: " << nclusters
        << std::endl;
    }

    bool first = true;
    for( const auto& [key,cluster]:range_adaptor(cluster_range))
    {

      Cluster cluster_struct;
      cluster_struct.layer = layer;
      cluster_struct.tile = tile;

      ++m_container->n_detector_clusters[detid];

      // per view clusters
      if( detid < 8 ) ++m_container->n_phi_clusters;
      else ++m_container->n_z_clusters;

      // find associated hits
      const auto range = m_cluster_hit_map->getHits(key);
      cluster_struct.size = std::distance( range.first, range.second );

      // loop over assiciated hits
      unsigned int adc_max = 0;
      for( auto iter = range.first; iter != range.second; ++iter )
      {
        const auto& [cluskey,hitkey] = *iter;
        assert( cluskey == key );

        // get strip
        const auto strip = MicromegasDefs::getStrip( hitkey );

        // get associated hit
        const auto hit = hitset->getHit( hitkey );
        assert( hit );

        const auto adc = hit->getAdc();
        if( adc > adc_max )
        {
          adc_max = adc;
          cluster_struct.strip = strip;
        }

        // get adc, remove pedestal, increment total charge
        const double pedestal = m_use_default_pedestal ?
          m_default_pedestal:
          m_calibration_data.get_pedestal_mapped(hitsetkey, strip);
        cluster_struct.charge += (adc-pedestal);

      }

      if( first )
      {
        m_container->first_cluster_strip[detid] = cluster_struct.strip;
        first = false;
      }

      if( min_cluster_size == 0 || cluster_struct.size < min_cluster_size ) { min_cluster_size = cluster_struct.size; }
      if( min_cluster_charge == 0 || cluster_struct.charge < min_cluster_charge ) { min_cluster_charge = cluster_struct.charge; }

      // cluster region
      cluster_struct.region = cluster_struct.strip/64;

      // increment number of clusters in region (0 to 63)
      int region_id = cluster_struct.region+4*detid;

      ++m_container->n_region_clusters[region_id];

      // local position
      cluster_struct.x_local = cluster->getLocalX();
      cluster_struct.y_local = cluster->getLocalY();

      // global position
      auto globalPosition = m_tGeometry->getGlobalPosition(key, cluster);
      cluster_struct.x = globalPosition.x();
      cluster_struct.y = globalPosition.y();
      cluster_struct.z = globalPosition.z();


      m_container->clusters.push_back(cluster_struct);

    }

    // store min cluster size
    m_container->min_cluster_size[detid] = min_cluster_size;

    // store min cluster charge
    m_container->min_cluster_charge[detid] = min_cluster_charge;

  }

}
