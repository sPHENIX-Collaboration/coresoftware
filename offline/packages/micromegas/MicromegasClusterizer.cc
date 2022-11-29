/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasClusterizer.h"
#include "MicromegasDefs.h"
#include "CylinderGeomMicromegas.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>           // for PHG4CylinderGeom

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainerv4.h>        // for TrkrCluster
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssocv3.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                     // for PHIODataNode
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>                   // for PHNodeIterator
#include <phool/PHObject.h>                         // for PHObject

#include <Eigen/Dense>

#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <cstdint>                                 // for uint16_t
#include <iterator>                                 // for distance
#include <map>                                      // for _Rb_tree_const_it...
#include <utility>                                  // for pair, make_pair
#include <vector>


namespace
{
  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  // streamers
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const Acts::Vector3& vector )
  {
    out << "( " << vector[0] << "," << vector[1] << "," << vector[2] << ")";
    return out;
  }

  // streamers
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const Acts::Vector2& vector )
  {
    out << "( " << vector[0] << "," << vector[1] << ")";
    return out;
  }

  // streamers
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const TVector3& vector )
  {
    out << "( " << vector.x() << "," << vector.y() << "," << vector.z() << ")";
    return out;
  }

  // streamers
  [[maybe_unused]] inline std::ostream& operator << (std::ostream& out, const TVector2& vector )
  {
    out << "( " << vector.X() << "," << vector.Y() << ")";
    return out;
  }

}

//_______________________________________________________________________________
MicromegasClusterizer::MicromegasClusterizer(const std::string &name, const std::string& detector)
  : SubsysReco(name)
  , m_detector( detector )
{}

//_______________________________________________________________________________
int MicromegasClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  assert( dstNode );

  // Create the Cluster node if missing
  auto trkrClusterContainer = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrClusterContainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto trkrNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if(!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    trkrClusterContainer = new TrkrClusterContainerv4;
    auto TrkrClusterContainerNode = new PHIODataNode<PHObject>(trkrClusterContainer, "TRKR_CLUSTER", "PHObject");
    trkrNode->addNode(TrkrClusterContainerNode);
  }

  // create cluster to hit association node, if missing
  auto trkrClusterHitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
  if(!trkrClusterHitAssoc)
  {
    PHNodeIterator dstiter(dstNode);
    auto trkrNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if(!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    trkrClusterHitAssoc = new TrkrClusterHitAssocv3;
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(trkrClusterHitAssoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    trkrNode->addNode(newNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________________________________________________
int MicromegasClusterizer::process_event(PHCompositeNode *topNode)
{

  // geometry
  PHG4CylinderGeomContainer* geonode = nullptr;
  for( const std::string& geonodename: {"CYLINDERGEOM_" + m_detector + "_FULL", "CYLINDERGEOM_" + m_detector } )
  { if(( geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str()) )) break; }
  assert(geonode);

  // hitset container
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert( trkrhitsetcontainer );

  // cluster container
  auto trkrClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert( trkrClusterContainer );

  // cluster-hit association
  auto trkrClusterHitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  assert( trkrClusterHitAssoc );

  // geometry
  auto acts_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( acts_geometry );

  // loop over micromegas hitsets
  const auto hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId);
  for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
  {
    // get hitset, key and layer
    TrkrHitSet* hitset = hitset_it->second;
    const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;
    const auto layer = TrkrDefs::getLayer(hitsetkey);
    const auto tileid = MicromegasDefs::getTileId(hitsetkey);

    // get micromegas geometry object
    const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(geonode->GetLayerGeom(layer));
    assert(layergeom);

    // get micromegas acts surface
    const auto acts_surface = acts_geometry->maps().getMMSurface( hitsetkey);
    if( !acts_surface )
    {
      std::cout
        << "MicromegasClusterizer::process_event -"
        << " could not find surface for layer " << (int) layer << " tile: " << (int) tileid
        << " skipping hitset"
        << std::endl;
      continue;
    }

    /*
     * get segmentation type, layer thickness, strip length and pitch.
     * They are used to calculate cluster errors
     */
    const auto segmentation_type = layergeom->get_segmentation_type();
    const double pitch = layergeom->get_pitch();
    const double strip_length = layergeom->get_strip_length( tileid, acts_geometry );
    
    // keep a list of ranges corresponding to each cluster
    using range_list_t = std::vector<TrkrHitSet::ConstRange>;
    range_list_t ranges;

    // loop over hits
    const auto hit_range = hitset->getHits();

    // keep track of first iterator of runing cluster
    auto begin = hit_range.first;

    // keep track of previous strip
    uint16_t previous_strip = 0;
    bool first = true;

    for( auto hit_it = hit_range.first; hit_it != hit_range.second; ++hit_it )
    {

      // get hit key
      const auto hitkey = hit_it->first;

      // get strip number
      const auto strip = MicromegasDefs::getStrip( hitkey );

      if( first )
      {

        previous_strip = strip;
        first = false;
        continue;

      } else if( strip - previous_strip > 1 ) {

        // store current cluster range
        ranges.push_back( std::make_pair( begin, hit_it ) );

        // reinitialize begin of next cluster range
        begin = hit_it;

      }

      // update previous strip
      previous_strip = strip;

    }

    // store last cluster
    if( begin != hit_range.second ) ranges.push_back( std::make_pair( begin, hit_range.second ) );

    // initialize cluster count
    int cluster_count = 0;
    int strip_count = 0;
    // loop over found hit ranges and create clusters
    for( const auto& range : ranges )
    {
      strip_count++;
      // create cluster key and corresponding cluster
      const auto ckey = TrkrDefs::genClusKey( hitsetkey, cluster_count++ );

      TVector2 local_coordinates;
      double weight_sum = 0;

      // needed for proper error calculation
      // it is either the sum over z, or phi, depending on segmentation
      double coord_sum = 0;
      double coordsquare_sum = 0;

      // also store adc value
      unsigned int adc_sum = 0;
      strip_count = 0;
      // loop over constituting hits
      for( auto hit_it = range.first; hit_it != range.second; ++hit_it )
      {
	strip_count++;
        // get hit key
        const auto hitkey = hit_it->first;
        const auto hit = hit_it->second;

        // associate cluster key to hit key
        trkrClusterHitAssoc->addAssoc(ckey, hitkey );

        // get strip number
        const auto strip = MicromegasDefs::getStrip( hitkey );

        // get adc, remove pedestal
        /* pedestal should be the same as the one used in PHG4MicromegasDigitizer */
        static constexpr double pedestal = 74.6;
        const double weight = double(hit->getAdc()) - pedestal;

        // increment cluster adc
        adc_sum += hit->getAdc();

        // get strip local coordinate and update relevant sums
        const auto strip_local_coordinate = layergeom->get_local_coordinates( tileid, acts_geometry, strip );
        local_coordinates += strip_local_coordinate*weight;
        switch( segmentation_type )
        {
          case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
          {

            coord_sum += strip_local_coordinate.X()*weight;
            coordsquare_sum += square(strip_local_coordinate.X())*weight;
            break;
          }

          case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
          {
            coord_sum += strip_local_coordinate.Y()*weight;
            coordsquare_sum += square(strip_local_coordinate.Y())*weight;
            break;
          }
        }

        weight_sum += weight;

      }

      local_coordinates *= (1./weight_sum);

      // dimension and error in r, rphi and z coordinates
      static const float invsqrt12 = 1./std::sqrt(12);
      static constexpr float error_scale_phi = 1.6;
      static constexpr float error_scale_z = 0.8;

      auto coord_cov = coordsquare_sum/weight_sum - square( coord_sum/weight_sum );
      auto coord_error_sq = coord_cov/weight_sum;

      // local errors (x is along rphi, y is along z)
      double error_sq_x = 0;
      double error_sq_y = 0;
      switch( segmentation_type )
      {
        case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
        {
          if( coord_error_sq == 0 ) coord_error_sq = square(pitch)/12;
          else coord_error_sq *= square(error_scale_phi);
          error_sq_x = coord_error_sq;
          error_sq_y = square(strip_length*invsqrt12);
          break;
        }
        
        case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
        {
          if( coord_error_sq == 0 ) coord_error_sq = square(pitch)/12;
          else coord_error_sq *= square(error_scale_z);
          error_sq_x = square(strip_length*invsqrt12);
          error_sq_y = coord_error_sq;
          break;
        }
      }

      
      if(m_cluster_version==3)
      {
        auto cluster = std::make_unique<TrkrClusterv3>();
        cluster->setAdc( adc_sum );
        cluster->setLocalX(local_coordinates.X());
        cluster->setLocalY(local_coordinates.Y());
        
        // assign errors
        cluster->setActsLocalError(0,0, error_sq_x);
        cluster->setActsLocalError(1,1, error_sq_y);
        cluster->setActsLocalError(0,1, 0);
        cluster->setActsLocalError(1,0, 0);

        // add to container
        trkrClusterContainer->addClusterSpecifyKey( ckey, cluster.release() );

      } else if(m_cluster_version==4) {

        auto cluster = std::make_unique<TrkrClusterv4>();
        cluster->setAdc( adc_sum );
        cluster->setLocalX(local_coordinates.X());
        cluster->setLocalY(local_coordinates.Y());

        // store cluster size
        switch( segmentation_type )
        {
          case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
          {
            cluster->setPhiSize(strip_count);
            cluster->setZSize(1);
            break;
          }

          case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
          {
            cluster->setPhiSize(1);
            cluster->setZSize(strip_count);
            break;
          }
        }
        // add to container
        trkrClusterContainer->addClusterSpecifyKey( ckey, cluster.release() );
      }

    }

  }
  // done
  return Fun4AllReturnCodes::EVENT_OK;
}
