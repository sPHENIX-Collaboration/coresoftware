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
#include <phool/PHNode.h>  
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TVector3.h>

#include <cassert>
#include <numeric>

namespace
{
  // convenient square function
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  // bind angle to [-M_PI,+M_PI[. This is useful to avoid edge effects when making the difference between two angles
  template<class T>
    inline T bind_angle( const T& angle )
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

      // distribute charge among adjacent strips
      // todo: will need to split charge for different segments along drift space
      // calculate the fired strips for each segment, including diffusion
      // get the complete cluster by resumming
      int tileid;
      charge_list_t fractions;
      std::tie(tileid, fractions) = distribute_charge( layergeom, world_mid, layergeom->get_pitch()/5 );
      if( tileid < 0 || fractions.empty() ) continue;

      // make sure fractions adds up to unity
      if( Verbosity() > 0 )
      {
        const auto sum = std::accumulate( fractions.begin(), fractions.end(), double( 0 ), 
          []( double value, const charge_pair_t& pair ) { return value + pair.second; } );
        std::cout << "PHG4MicromegasHitReco::process_event - sum: " << sum << std::endl;
      }

      // create hitset
      TrkrDefs::hitsetkey hitsetkey = MicromegasDefs::genHitSetKey( layer, tileid );
      auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);
      
      // generate the key for this hit
      // loop over strips in list
      for( const auto pair:fractions )
      {
        
        // get strip and bound check
        const int strip = pair.first;        
        if( strip < 0 || strip >= layergeom->get_strip_count( tileid ) ) continue;
        
        // get hit from hitset
        TrkrDefs::hitkey hitkey = MicromegasDefs::genHitKey(strip);
        TrkrHit* hit = hitset_it->second->getHit(hitkey);
        if( !hit )
        {
          // create hit and insert in hitset
          hit = new TrkrHit;
          hitset_it->second->addHitSpecificKey(hitkey, hit);
        }

        // add energy from g4hit       
        const double fraction = pair.second;
        hit->addEnergy( g4hit->get_eion()*fraction );
        
        // associate this hitset and hit to the geant4 hit key
        hittruthassoc->addAssoc(hitsetkey, hitkey, g4hit_it->first);
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
    cylinder->set_pitch( is_first ? 25./256 : 50./256 );
    cylinder->identify( std::cout );
  }
}

//___________________________________________________________________________
PHG4MicromegasHitReco::charge_info_t PHG4MicromegasHitReco::distribute_charge( CylinderGeomMicromegas* layergeom, const TVector3& location, float sigma ) const
{
  
  // find tile and strip matching center position
  int tileid, stripnum;
  std::tie(tileid, stripnum) = layergeom->find_strip( location );

  // check tile and strip
  if( tileid < 0 || stripnum < 0 ) return std::make_pair( -1, charge_list_t() );

  // store pitch and radius
  const auto pitch = layergeom->get_pitch();
  const auto radius = layergeom->get_radius();
  
  // find relevant strip indices
  const auto strip_count = layergeom->get_strip_count( tileid );
  const auto stripnum_min = std::clamp<int>( stripnum - 5.*sigma/pitch - 1, 0, strip_count );
  const auto stripnum_max = std::clamp<int>( stripnum + 5*sigma/pitch + 1, 0, strip_count );
  
  // prepare charge list
  charge_list_t charge_list;

  const auto phi = std::atan2( location.y(), location.x() );
  
  // loop over strips
  static const double sqrt2 = std::sqrt(2.);
  for( int strip = stripnum_min; strip <= stripnum_max; ++strip )
  {    
    // get strip center
    const TVector3 strip_location = layergeom->get_world_coordinate( tileid, strip );
    const auto phi_strip = std::atan2( strip_location.y(), strip_location.x() );
  
    // find relevant strip coordinate with respect to location
    const auto x_loc = layergeom->get_segmentation_type() == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ? 
      radius * bind_angle( phi_strip - phi ):
      strip_location.z() - location.z();
        
    // calculate charge fraction
//     // this corresponds to integrating the gaussian of width sigma from x_loc - pitch/2 to x_loc+pitch/2
//     const float fraction = (std::erf( (x_loc + pitch/2)/(sqrt2*sigma) ) - std::erf( (x_loc - pitch/2)/(sqrt2*sigma) ))/2;
    
    // this corresponds to zigzag strips with full overlap between one strip and its neighbors
    /* todo: 
     * - double check formulas
     * - match actual zigzag design from micromegas
     * - compare to GEM implementation
     */
    auto gaus = []( const double x, const double sigma ) { return std::exp( -square(x/sigma)/2 )/(sigma*std::sqrt(2*M_PI)); };
    const float fraction = 
      (pitch - x_loc)*(std::erf(x_loc/(sqrt2*sigma)) - erf((x_loc-pitch)/(sqrt2*sigma)))/(pitch*2)
      + (pitch + x_loc)*(erf((x_loc+pitch)/(sqrt2*sigma)) - std::erf(x_loc/(sqrt2*sigma)))/(pitch*2) 
      + (gaus(x_loc-pitch, sigma) - gaus(x_loc, sigma))*square(sigma)/pitch
      + (gaus(x_loc+pitch, sigma) - gaus(x_loc, sigma))*square(sigma)/pitch; 
          
    // store
    charge_list.push_back( std::make_pair( strip, fraction ) );
  }
    
  return std::make_pair( tileid, charge_list ); 
}
