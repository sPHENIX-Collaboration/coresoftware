/*!
 * \file PHG4MicromegasHitReco.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasHitReco.h"

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TVector3.h>

#include <cassert>

namespace
{
  // convenient square function
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  // bind angle to [-M_PI,+M_PI[. This is useful to avoid edge effects when making the difference between two angles
  template<class T>
    inline constexpr T bind_angle( const T& angle )
  {
    if( angle >= M_PI ) return angle - 2*M_PI;
    else if( angle < -M_PI ) return angle + 2*M_PI;
    else return angle;
  }

}

//___________________________________________________________________________
PHG4MicromegasHitReco::PHG4MicromegasHitReco(const std::string &name, const std::string& detector)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , m_detector(detector)
{ SetDefaultParameters(); }

//___________________________________________________________________________
int PHG4MicromegasHitReco::InitRun(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);

  // setup tiles
  setup_tiles( topNode );

  // get dst node
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << PHWHERE << "creating TRKR_HITSET." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    // create container and add to the tree
    hitsetcontainer = new TrkrHitSetContainer();
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }

  // create hit truth association if needed
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    std::cout << PHWHERE << "creating TRKR_HITTRUTHASSOC." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    hittruthassoc = new TrkrHitTruthAssoc();
    auto newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
    trkrnode->addNode(newNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
int PHG4MicromegasHitReco::process_event(PHCompositeNode *topNode)
{
  // load relevant nodes
  // G4Hits
  const std::string g4hitnodename = "G4HIT_" + m_detector;
  auto g4hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, g4hitnodename);
  assert(g4hitcontainer);

  // geometry
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  assert(geonode);

  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  // Get the TrkrHitTruthAssoc node
  auto hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  assert(hittruthassoc);

  // loop over layers in the g4hit container
  auto layer_range = g4hitcontainer->getLayers();
  for( auto layer_it = layer_range.first; layer_it != layer_range.second; ++layer_it )
  {
    // get layer
    const auto layer = *layer_it;

    // get relevant geometry
    auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(geonode->GetLayerGeom(layer));
    assert( layergeom );

    // get corresponding hits
    PHG4HitContainer::ConstRange g4hit_range = g4hitcontainer->getHits(layer);

    // loop over hits
    for( auto g4hit_it = g4hit_range.first; g4hit_it != g4hit_range.second; ++g4hit_it )
    {

      // get hit
      PHG4Hit* g4hit = g4hit_it->second;

      // check time window
      if(g4hit->get_t(0) > m_tmax) continue;
      if(g4hit->get_t(1) < m_tmin) continue;

      // get world coordinates
      TVector3 world_in( g4hit->get_x(0), g4hit->get_y(0), g4hit->get_z(0) );
      TVector3 world_out( g4hit->get_x(1), g4hit->get_y(1), g4hit->get_z(1) );
      auto world_mid = (world_in + world_out)*0.5;

      // get detector id and strip number from geometry
      const auto [tileid, stripnum] = layergeom->find_strip( world_mid );

      // check tile and strip
      if( tileid < 0 ) continue;
      if( stripnum < 0 ) continue;

      // for now, we create one hit per g4hit
      // TODO: implement proper smearing, lorentz angle, charge sharing among adjacent zig-zag strips, etc.
      // generate hitset key for this tile, and get corresponding hitset
      TrkrDefs::hitsetkey hitsetkey = MicromegasDefs::genHitSetKey( layer, tileid );
      auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);

      // generate the key for this hit
      TrkrDefs::hitkey hitkey = MicromegasDefs::genHitKey(stripnum);

      // get hit from hitset
      TrkrHit* hit = hitset_it->second->getHit(hitkey);
      if( !hit )
      {
        // create hit and insert in hitset
        hit = new TrkrHit;
        hitset_it->second->addHitSpecificKey(hitkey, hit);
      }

      // add energy from g4hit
      hit->addEnergy( g4hit->get_eion() );

      // associate this hitset and hit to the geant4 hit key
      hittruthassoc->addAssoc(hitsetkey, hitkey, g4hit_it->first);

      // debugging
      if(Verbosity() > 0)
      {
        auto hit_world_coord = layergeom->get_world_coordinate( tileid, stripnum );
        switch( layergeom->get_segmentation_type() )
        {

          case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
          {
            auto truth_phi = std::atan2( world_mid.y(), world_mid.x() );
            auto hit_phi = std::atan2( hit_world_coord.y(), hit_world_coord.x() );
            std::cout << "PHG4MicromegasHitReco::process_event - layer: " << layer << " rphi difference: " << bind_angle( hit_phi - truth_phi )*layergeom->get_radius() << std::endl;
            break;
          }

          case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
          {
            auto truth_z = world_mid.z();
            auto hit_z = hit_world_coord.z();
            std::cout << "PHG4MicromegasHitReco::process_event - layer: " << layer << "    z difference: " << hit_z - truth_z << std::endl;
            break;
          }

        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void PHG4MicromegasHitReco::SetDefaultParameters()
{
  // default timing window
  // TODO: allow to modify this via PHParameterInterface
  m_tmin = -5000;
  m_tmax = 5000;
}

//___________________________________________________________________________
void PHG4MicromegasHitReco::setup_tiles(PHCompositeNode* topNode)
{

  // get geometry
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geonode = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename);
  if (!geonode)
  {
    std::cout << PHWHERE << "Could not locate geometry node " << geonodename << std::endl;
    exit(1);
  }

  // get cylinders
  PHG4CylinderGeomContainer::ConstRange range = geonode->get_begin_end();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    std::cout << "PHG4MicromegasHitReco::setup_tiles - processing layer " << iter->first << std::endl;
    auto cylinder = static_cast<CylinderGeomMicromegas*>(iter->second);

    // assign tiles
    cylinder->set_tiles( m_tiles );

    // asign segmentation type and pitch
    // assume first layer in phi, other(s) are z
    const bool is_first( iter == range.first );
    cylinder->set_segmentation_type( is_first ?
      MicromegasDefs::SegmentationType::SEGMENTATION_PHI :
      MicromegasDefs::SegmentationType::SEGMENTATION_Z );

    // pitch
    /* they correspond to 256 channels along the phi direction, and 256 along the z direction, assuming 25x50 tiles */
    cylinder->set_pitch( is_first ? 0.0977 : 0.195 );

  }
}

