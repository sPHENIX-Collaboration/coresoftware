/*!
 * \file MicromegasClusterizer.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasClusterizer.h"
#include "MicromegasDefs.h"
#include "CylinderGeomMicromegas.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <Eigen/Dense>

#include <TVector3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>

namespace
{
  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
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

    trkrClusterContainer = new TrkrClusterContainer();
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

    trkrClusterHitAssoc = new TrkrClusterHitAssoc();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(trkrClusterHitAssoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    trkrNode->addNode(newNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________________________________________________
int MicromegasClusterizer::process_event(PHCompositeNode *topNode)
{

  // geometry
  const std::string geonodename = "CYLINDERGEOM_" + m_detector;
  auto geonode =  findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
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

  // loop over micromegas hitsets
  const auto hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId);
  for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
  {

    // get hitset, key and layer
    TrkrHitSet* hitset = hitset_it->second;
    const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;
    const auto layer = TrkrDefs::getLayer(hitsetkey);
    const auto tileid = MicromegasDefs::getTileId(hitsetkey);

    // get geometry object
    const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(geonode->GetLayerGeom(layer));
    assert(layergeom);

    /*
     * get segmentation type, layer thickness, strip length and pitch.
     * They are used to calculate cluster errors
     */
    const auto segmentation_type = layergeom->get_segmentation_type();
    const double thickness = layergeom->get_thickness();
    const double pitch = layergeom->get_pitch();
    const double strip_length = layergeom->get_strip_length( tileid );

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

    // loop over found hit ranges and create clusters
    for( const auto& range : ranges )
    {

      // create cluster key and corresponding cluster
      const auto cluster_key = MicromegasDefs::genClusterKey( hitsetkey, cluster_count++ );
      auto cluster = (trkrClusterContainer->findOrAddCluster(cluster_key))->second;

      TVector3 world_coordinates;
      double adc_sum = 0;

      // loop over constituting hits
      for( auto hit_it = range.first; hit_it != range.second; ++hit_it )
      {

        // get hit key
        const auto hitkey = hit_it->first;
        const auto hit = hit_it->second;

        // associate cluster key to hit key
        trkrClusterHitAssoc->addAssoc(cluster_key, hitkey );

        // get strip number
        const auto strip = MicromegasDefs::getStrip( hitkey );

        // get adc, remove pedestal
        /* pedestal should be the same as the one used in PHG4MicromegasDigitizer */
        static constexpr double pedestal = 74.6;
        const double weight = double(hit->getAdc()) - pedestal;

        // get strip world coordinate
        world_coordinates += layergeom->get_world_coordinate( tileid, strip )*weight;
        adc_sum += weight;

      }

      // cluster position
      cluster->setPosition( 0, world_coordinates.x()/adc_sum );
      cluster->setPosition( 1, world_coordinates.y()/adc_sum );
      cluster->setPosition( 2, world_coordinates.z()/adc_sum );
      cluster->setGlobal();

      // dimension and error in r, rphi and z coordinates
      static const float invsqrt12 = 1./std::sqrt(12);

      using matrix_t = Eigen::Matrix<float, 3, 3>;
      matrix_t dimension = matrix_t::Zero();
      matrix_t error = matrix_t::Zero();

      const auto size = std::distance( range.first, range.second );
      switch( segmentation_type )
      {
        case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
        dimension(0,0) = square(0.5*thickness);
        dimension(1,1) = square( 0.5*pitch*size );
        dimension(2,2) = square( 0.5*strip_length );

        error(0,0) = square(thickness*invsqrt12);
        error(1,1) = square( pitch*invsqrt12 );
        error(2,2) = square( strip_length*invsqrt12 );
        break;

        case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
        dimension(0,0) = square(0.5*thickness);
        dimension(1,1) = square( 0.5*strip_length );
        dimension(2,2) = square( 0.5*pitch*size );

        error(0,0) = square(thickness*invsqrt12);
        error(1,1) = square( strip_length*invsqrt12 );
        error(2,2) = square( pitch*invsqrt12 );
        break;
      }

      // rotate and save
      matrix_t rotation = matrix_t::Identity();
      const double phi = std::atan2( world_coordinates.y(), world_coordinates.x() );
      const double cosphi = std::cos(phi);
      const double sinphi = std::sin(phi);
      rotation(0,0) = cosphi;
      rotation(0,1) = -sinphi;
      rotation(1,0) = sinphi;
      rotation(1,1) = cosphi;

      // rotate dimension and error
      dimension = rotation*dimension*rotation.transpose();
      error = rotation*error*rotation.transpose();

      // assign to cluster
      for( int i = 0; i<3; ++i )
        for( int j = 0; j<3; ++j )
      {
        cluster->setSize( i, j, dimension(i,j) );
        cluster->setError( i, j, error(i,j) );
      }

    }

  }

  // done
  return Fun4AllReturnCodes::EVENT_OK;
}

