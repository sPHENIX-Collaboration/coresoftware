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

  // do clusterization
  // for now each found hit is converted into a cluster
  // TODO: implement real cluster

  // loop over micromegas hitsets
  const auto hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::micromegasId);
  for( auto hitset_it = hitset_range.first; hitset_it != hitset_range.second; ++hitset_it )
  {

    // get hitset, key and layer
    TrkrHitSet* hitset = hitset_it->second;
    const TrkrDefs::hitsetkey hitsetkey = hitset_it->first;
    const auto layer = TrkrDefs::getLayer(hitsetkey);
    const auto tile = MicromegasDefs::getTileId(hitsetkey);

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
    const double strip_length = layergeom->get_strip_length( tile );

    // initialize cluster count
    int cluster_count = 0;

    // loop over hits
    const auto hit_range = hitset->getHits();
    for( auto hit_it = hit_range.first; hit_it != hit_range.second; ++hit_it )
    {

      // get hit key
      const auto hitkey = hit_it->first;

      // get strip number
      const auto strip = MicromegasDefs::getStrip( hitkey );

      // get strip world coordinate
      TVector3 world_coordinates = layergeom->get_world_coordinate( tile, strip );

      // create cluster key
      const auto cluster_key = MicromegasDefs::genClusterKey( hitsetkey, cluster_count++ );

      // associate cluster key to hit key
      trkrClusterHitAssoc->addAssoc(cluster_key, hitkey );

      // create new cluster of this key
      auto cluster = (trkrClusterContainer->findOrAddCluster(cluster_key))->second;

      // cluster position
      cluster->setPosition( 0, world_coordinates.x() );
      cluster->setPosition( 1, world_coordinates.y() );
      cluster->setPosition( 2, world_coordinates.z() );
      cluster->setGlobal();

      // dimension and error in r, rphi and z coordinates
      static const float invsqrt12 = 1./std::sqrt(12);

      using matrix_t = Eigen::Matrix<float, 3, 3>;
      matrix_t dimension = matrix_t::Zero();
      matrix_t error = matrix_t::Zero();

      switch( segmentation_type )
      {
        case MicromegasDefs::SegmentationType::SEGMENTATION_PHI:
        dimension(0,0) = square(0.5*thickness);
        dimension(1,1) = square( 0.5*pitch );
        dimension(2,2) = square( 0.5*strip_length );

        error(0,0) = square(thickness*invsqrt12);
        error(1,1) = square( pitch*invsqrt12 );
        error(2,2) = square( strip_length*invsqrt12 );
        break;

        case MicromegasDefs::SegmentationType::SEGMENTATION_Z:
        dimension(0,0) = square(0.5*thickness);
        dimension(1,1) = square( 0.5*strip_length );
        dimension(2,2) = square( 0.5*pitch );

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

