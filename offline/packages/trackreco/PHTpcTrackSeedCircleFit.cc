#include "PHTpcTrackSeedCircleFit.h"

#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/ActsTransformations.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <TF1.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, std::sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>

namespace
{
   
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // radius
  template<class T> T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

}

//____________________________________________________________________________..
PHTpcTrackSeedCircleFit::PHTpcTrackSeedCircleFit(const std::string &name):
 SubsysReco(name)
{}

//____________________________________________________________________________..
int PHTpcTrackSeedCircleFit::InitRun(PHCompositeNode *topnode)
{ 
  // get relevant nodes
  int ret = GetNodes( topnode );
  if( ret != Fun4AllReturnCodes::EVENT_OK ) return ret;

  return Fun4AllReturnCodes::EVENT_OK; }

//____________________________________________________________________________..
int PHTpcTrackSeedCircleFit::process_event(PHCompositeNode*)
{
  
  // _track_map contains the TPC seed track stubs
  // We want to associate these TPC track seeds with a collision vertex
  // All we need is to project the TPC clusters in Z to the beam line.
  // The TPC track seeds are given a position that is the PCA of the line and circle fit to the beam line

  if(Verbosity() > 0)
    std::cout << PHWHERE << " TPC track map size " << _track_map->size()  << std::endl;

  unsigned int track_key = 0;
  // loop over the TPC track seeds
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      TrackSeed* tracklet_tpc = *phtrk_iter;
      
      if (Verbosity() > 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << track_key
	    << ": nhits: " << tracklet_tpc-> size_cluster_keys()
	    << ": pT: " << tracklet_tpc->get_pt()
	    << ": phi: " << tracklet_tpc->get_phi(_cluster_map, _surfmaps, _tGeometry)
	    << ": eta: " << tracklet_tpc->get_eta()
	    << std::endl;
	}

      if(tracklet_tpc->size_cluster_keys() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough TPC clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}
      
      /// Start at layer 7
      tracklet_tpc->lineFit(_cluster_map, _surfmaps, _tGeometry, 7, 58);
      tracklet_tpc->circleFitByTaubin(_cluster_map, _surfmaps, _tGeometry, 7, 58);
   
      if(Verbosity() > 5)
      {
        std::cout << " new mom " <<  tracklet_tpc->get_p() <<  "  new eta " <<  tracklet_tpc->get_eta()
          << " new phi " << tracklet_tpc->get_phi(_cluster_map, _surfmaps, _tGeometry) * 180.0 / M_PI << std::endl;
      }

      track_key++;
    }  // end loop over TPC track seeds

    if(Verbosity() > 0)
      std::cout << " Final track map size " << _track_map->size() << std::endl;

    if (Verbosity() > 0)
      std::cout << "PHTpcTrackSeedCircleFit::process_event(PHCompositeNode *topNode) Leaving process_event" << std::endl;

    return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcTrackSeedCircleFit::End(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }

int  PHTpcTrackSeedCircleFit::GetNodes(PHCompositeNode* topNode)
{
  _surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  if(!_surfmaps)
    {
      std::cout << PHWHERE << "Error, can't find acts surface maps" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if(_use_truth_clusters)
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  else
    {
      _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    }
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  if( _dcc )
  { std::cout << "PHTpcTrackSeedCircleFit::get_Nodes  - found static TPC distortion correction container" << std::endl; }

  _track_map = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

  
Acts::Vector3 PHTpcTrackSeedCircleFit::getGlobalPosition( TrkrDefs::cluskey key, TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  auto globalpos = _surfmaps->getGlobalPosition(key, cluster, _tGeometry);

    // ADF: in streaming mode we need to add a step here to take care of the fact that we do not know the crossing yet
  // possibly move the track to point at z=0 to make distortion corrections (circularize the track) then move it back after the fit?

  // check if TPC distortion correction are in place and apply
  if(_dcc) { globalpos = _distortionCorrection.get_corrected_position( globalpos, _dcc ); }

  return globalpos;
}

