#include "MicromegasEvaluator_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/PHRandomSeed.h>
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

  /// square
  template<class T> inline constexpr T square( T x ) { return x*x; }

  //* converninece trait for underlying type
  template<class T>
    using underlying_type_t = typename std::underlying_type<T>::type;

  //* convert an strong type enum to integral type
  template<class T>
    constexpr underlying_type_t<T>
    to_underlying_type(T value) noexcept
  { return static_cast<underlying_type_t<T>>(value);}

  /// create cluster struct from svx cluster
  MicromegasEvaluator_hp::HitStruct create_hit( TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey, TrkrHit* hit )
  {
    MicromegasEvaluator_hp::HitStruct hit_struct;
    hit_struct._layer = TrkrDefs::getLayer(hitsetkey);
    hit_struct._tile = MicromegasDefs::getTileId(hitsetkey);
    hit_struct._strip = MicromegasDefs::getStrip(hitkey);
    hit_struct._energy = hit->getEnergy();
    hit_struct._adc = hit->getAdc();
    return hit_struct;
  }

}

//_____________________________________________________________________
void MicromegasEvaluator_hp::Container::Reset()
{ _hits.clear(); }

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
  const std::string geonodename = "CYLINDERGEOM_MICROMEGAS";
  const auto geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
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

  evaluate_g4hits();
  evaluate_hits();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int MicromegasEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // geometry
  const std::string geonodename = "CYLINDERGEOM_MICROMEGAS";
  m_geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  assert(m_geonode);

  // g4hit container
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");
  assert(m_g4hits_micromegas);

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(m_hitsetcontainer);

  // local container
  m_container = findNode::getClass<Container>(topNode, "MicromegasEvaluator_hp::Container");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void MicromegasEvaluator_hp::evaluate_g4hits()
{

  if( !( m_g4hits_micromegas && m_container ) ) return;

  // loop over layers in the g4hit container
  auto layer_range = m_g4hits_micromegas->getLayers();
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
      PHG4Hit* g4hit = g4hit_it->second;

      // get world coordinates
      TVector3 world_in( g4hit->get_x(0), g4hit->get_y(0), g4hit->get_z(0) );
      TVector3 world_out( g4hit->get_x(1), g4hit->get_y(1), g4hit->get_z(1) );
      auto location = (world_in + world_out)*0.5;

      // get tile and strip, skip if invalid
      int tileid, stripnum;
      std::tie(tileid, stripnum) = layergeom->find_strip( location );

      // check tile and strip
      if( tileid < 0 || stripnum < 0 ) continue;

      // create G4Hit struct
      G4HitStruct g4hit_struct;
      g4hit_struct._layer = layer;
      g4hit_struct._tile = tileid;
      g4hit_struct._eion = g4hit->get_eion();

      // copied from PHG4MicromegasHitReco
      static constexpr double electrons_per_gev = 3.73252e+07;
      static constexpr double gain = 2000;

      // generate number of primary electrons
      g4hit_struct._nprimary = gsl_ran_poisson(m_rng.get(), g4hit_struct._eion*electrons_per_gev);

      // calculate total number of electrons and assign
      uint ntot = 0;
      for( uint i = 0; i < g4hit_struct._nprimary; ++i )
      { ntot += gsl_ran_exponential(m_rng.get(), gain); }

      g4hit_struct._nelectron = ntot;

      // store
      m_container->addG4Hit( g4hit_struct );

    }
  }

}

//_____________________________________________________________________
void MicromegasEvaluator_hp::evaluate_hits()
{

  if( !( m_hitsetcontainer && m_container ) ) return;

  // clear array
  m_container->clearHits();

  // loop over micromegas and tpc hitsets
  for( auto id:{TrkrDefs::TrkrId::micromegasId, TrkrDefs::TrkrId::tpcId} )
  {

    const auto hitset_range = m_hitsetcontainer->getHitSets(id);
    for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
    {

      // get hitset, key and layer
      TrkrHitSet* hitset = hitset_it->second;
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

      }
    }

  }

}
