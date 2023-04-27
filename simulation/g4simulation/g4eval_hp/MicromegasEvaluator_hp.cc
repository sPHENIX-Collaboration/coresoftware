#include "MicromegasEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasMapping.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/PHRandomSeed.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <gsl/gsl_randist.h>

//_____________________________________________________________________
namespace
{

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

  //! radius
  template<class T> inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  //! radius
  float get_r( PHG4Hit* hit, int i )
  {  return get_r( hit->get_x(i), hit->get_y(i) ); }

  //! create g4hit struct from G4Hit
  MicromegasEvaluator_hp::G4HitStruct create_g4hit( PHG4Hit* g4hit )
  {
    MicromegasEvaluator_hp::G4HitStruct g4hitstruct;
    g4hitstruct._edep = g4hit->get_edep();
    g4hitstruct._eion = g4hit->get_eion();
    g4hitstruct._t = g4hit->get_avg_t();
    g4hitstruct._x = g4hit->get_avg_x();
    g4hitstruct._y = g4hit->get_avg_y();
    g4hitstruct._z = g4hit->get_avg_z();
    g4hitstruct._delta_r = get_r( g4hit, 1 ) - get_r( g4hit, 0 );
    return g4hitstruct;
  }

  //! create cluster struct from svx cluster
  MicromegasEvaluator_hp::HitStruct create_hit( TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey, TrkrHit* hit )
  {
    MicromegasEvaluator_hp::HitStruct hit_struct;
    hit_struct._detid = TrkrDefs::getTrkrId(hitsetkey);
    hit_struct._layer = TrkrDefs::getLayer(hitsetkey);
    hit_struct._tile = MicromegasDefs::getTileId(hitsetkey);
    hit_struct._strip = MicromegasDefs::getStrip(hitkey);
    hit_struct._energy = hit->getEnergy();
    hit_struct._adc = hit->getAdc();
    return hit_struct;
  }
  
  // TVector3 streamer
  std::ostream& operator << (std::ostream& out, const TVector3& v )
  {
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }
  
}

//_____________________________________________________________________
void MicromegasEvaluator_hp::Container::Reset()
{
  _tiles.clear();
  _g4hits.clear();
  _hits.clear();
}

//_____________________________________________________________________
MicromegasEvaluator_hp::TileStruct& MicromegasEvaluator_hp::Container::findTile( uint layer, uint tile )
{
  auto iter = std::find_if( _tiles.begin(), _tiles.end(),
    [layer, tile]( const TileStruct& tilestruct )
  { return tilestruct._layer == layer && tilestruct._tile == tile; } );

  if( iter != _tiles.end() ) return *iter;
  else {
    _tiles.emplace_back( layer, tile );
    return _tiles.back();
  }

}

//_____________________________________________________________________
MicromegasEvaluator_hp::MicromegasEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{
  // initialize rng
  const uint seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );
}

//_____________________________________________________________________
int MicromegasEvaluator_hp::Init(PHCompositeNode* topNode )
{
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MicromegasEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "MicromegasEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  auto newNode = new PHIODataNode<PHObject>( new Container, "MicromegasEvaluator_hp::Container","PHObject");
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasEvaluator_hp::InitRun(PHCompositeNode* topNode )
{
  // load geometry
  PHG4CylinderGeomContainer* geonode = nullptr;
  for( const auto& geonodename: {"CYLINDERGEOM_MICROMEGAS_FULL", "CYLINDERGEOM_MICROMEGAS" } )
  { if(( geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename))) break; }
  assert(geonode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;
  if( m_container ) m_container->Reset();

  if( m_flags & EvalG4Hits ) evaluate_g4hits();
  if( m_flags & EvalHits ) evaluate_hits();

  static bool first = true;
  if( first )
  {
    first = false;
    if( m_flags & PrintGeometry ) print_micromegas_geometry(); 
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int MicromegasEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // geometry
  const std::string geonodename = "CYLINDERGEOM_MICROMEGAS_FULL";
  m_geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  assert(m_geonode);

  // g4hit container
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // local container
  m_container = findNode::getClass<Container>(topNode, "MicromegasEvaluator_hp::Container");
  assert(m_container);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void MicromegasEvaluator_hp::evaluate_g4hits()
{
  if( !( m_g4hits_micromegas && m_container ) )
  {
    std::cerr << "MicromegasEvaluator_hp::evaluate_g4hits - nodes not found" << std::endl;
    return;
  }

  // clear array
  m_container->clearG4Hits();

  // loop over layers in the g4hit container
  const auto layer_range = m_g4hits_micromegas->getLayers();

  TileStruct* current_tile = nullptr;
  for( auto layer_it = layer_range.first; layer_it != layer_range.second; ++layer_it )
  {

    // get layer
    const auto layer = *layer_it;

    // get relevant geometry
    auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_geonode->GetLayerGeom(layer));
    assert( layergeom );

    // get corresponding hits
    PHG4HitContainer::ConstRange g4hit_range = m_g4hits_micromegas->getHits(layer);

    // loop over hits
    for( auto g4hit_it = g4hit_range.first; g4hit_it != g4hit_range.second; ++g4hit_it )
    {
      // get hit
      auto g4hit = g4hit_it->second;

      const auto tileid = g4hit->get_property_int( PHG4Hit::prop_index_i );
      if( tileid < 0 ) continue;

      // get world coordinates
      TVector3 world_in( g4hit->get_x(0), g4hit->get_y(0), g4hit->get_z(0) );
      TVector3 world_out( g4hit->get_x(1), g4hit->get_y(1), g4hit->get_z(1) );

      const int stripnum = layergeom->find_strip_from_world_coords( tileid, m_tGeometry, (world_in + world_out)*0.5  );
      if( stripnum < 0 ) continue;

      // create G4Hit struct
      auto g4hitstruct = create_g4hit( g4hit );
      g4hitstruct._layer = layer;
      g4hitstruct._tile = tileid;

      // copied from PHG4MicromegasHitReco
      static constexpr double electrons_per_gev = 3.73252e+07;
      static constexpr double gain = 2000;

      // generate number of primary electrons
      g4hitstruct._nprimary = gsl_ran_poisson(m_rng.get(), g4hitstruct._eion*electrons_per_gev);

      // calculate total number of electrons and assign
      uint ntot = 0;
      for( uint i = 0; i < g4hitstruct._nprimary; ++i )
      { ntot += gsl_ran_exponential(m_rng.get(), gain); }

      g4hitstruct._nelectron = ntot;

      // store
      m_container->addG4Hit( g4hitstruct );

      // also fill tile information
      if( !( current_tile && current_tile->_layer == layer && int( current_tile->_tile ) == tileid ) )
      { current_tile = &m_container->findTile( layer, tileid ); }

      // fill
      current_tile->_edep_total += g4hitstruct._edep;
      current_tile->_eion_total += g4hitstruct._eion;

    }

  }

}

//_____________________________________________________________________
void MicromegasEvaluator_hp::evaluate_hits()
{

  if( !( m_hitsetcontainer && m_container ) ) return;

  // clear array
  m_container->clearHits();

  // keep track of current tile struct
  TileStruct* current_tile = nullptr;

  // loop over micromegas and tpc hitsets
  for( auto id:{TrkrDefs::TrkrId::micromegasId, TrkrDefs::TrkrId::tpcId} )
  // for( auto id:{TrkrDefs::TrkrId::micromegasId} )
  {

    const auto hitset_range = m_hitsetcontainer->getHitSets(id);
    for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
    {

      // get hitset, key and layer
      auto hitset = hitset_it->second;
      const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;

      // loop over hits
      const auto hit_range = hitset->getHits();

      for( auto hit_it = hit_range.first; hit_it != hit_range.second; ++hit_it )
      {

        // get key and hit
        const auto hitkey = hit_it->first;
        const auto hit = hit_it->second;

        // create hit struct
        auto hit_struct = create_hit( hitsetkey, hitkey, hit );

        // store
        m_container->addHit( hit_struct );

        // also fill tile information
        if( !( current_tile && current_tile->_layer == hit_struct._layer && current_tile->_tile == hit_struct._tile ) )
        { current_tile = &m_container->findTile( hit_struct._layer, hit_struct._tile ); }

        // fill
        current_tile->_electrons_total += hit_struct._energy;
        ++current_tile->_strips_total;

      }
    }

  }

}


//_____________________________________________________________________
void MicromegasEvaluator_hp::print_micromegas_geometry()
{
  std::cout << "MicromegasEvaluator_hp::print_micromegas_geometry" << std::endl;

  MicromegasMapping mapping;
  
  // loop over tiles
  for( uint tileid = 0; tileid < 8; ++tileid )
  {
    
    for( int layer:{55,56} )
    {
      // get relevant geometry
      auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_geonode->GetLayerGeom(layer));
      assert( layergeom );
      
      const auto thickness = layergeom->get_thickness();
      const auto strip_length = layergeom->get_strip_length( tileid, m_tGeometry);      
            
      for( const uint stripnum:{0,255} )
      {
        
        // get strip center, local coordinates
        const auto strip_center_local = layergeom->get_local_coordinates( tileid, m_tGeometry, stripnum );
        TVector3 strip_begin_world;
        TVector3 strip_end_world;
        
        // delta_z strip accounts for strip position along z with respect to middle of drift space
        const double delta_z_strip = thickness/2
          + 18e-4  // mesh
          + 120e-4 // amplification
          + 20e-4  // resist layer
          + 50e-4;  // kapton
        
        const auto segmentation = layergeom->get_segmentation_type();
        switch( segmentation )
        {
          case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
          {
            
            const TVector3 strip_begin( strip_center_local.X(), -strip_length/2, delta_z_strip );
            const TVector3 strip_end( strip_center_local.X(), strip_length/2, delta_z_strip );
            
            // convert to world coordinates and save
            strip_begin_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, strip_begin );
            strip_end_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, strip_end );            
            break;
          }
          
          case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
          {
            const TVector3 strip_begin( -strip_length/2, strip_center_local.Y(), -delta_z_strip );
            const TVector3 strip_end( strip_length/2, strip_center_local.Y(), -delta_z_strip );

            // convert to world coordinates and save
            strip_begin_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, strip_begin );
            strip_end_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, strip_end );            
            break;
          }
      
        }
       
        // get module name from layer and tileid
        const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, tileid );
        const auto module_name_sphenix = mapping.get_detname_sphenix_from_hitsetkey(hitsetkey);
        const auto module_name_saclay = mapping.get_detname_saclay_from_hitsetkey(hitsetkey);

        std::cout << layer << " " << tileid << " " << module_name_sphenix << " " << module_name_saclay << " " << stripnum << " " <<  strip_begin_world << std::endl;
        std::cout << layer << " " << tileid << " " << module_name_sphenix << " " << module_name_saclay << " " << stripnum << " " <<  strip_end_world << std::endl;
      }
      std::cout << std::endl;
    }
  }

}
