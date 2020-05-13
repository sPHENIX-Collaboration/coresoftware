/*!
 *  \file       PHTpcSeedFinder.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#include "PHTpcSeedFinder.h"
#include "PHTpcTrackerUtil.h"

#include <phool/PHLog.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>

float round(float var)
{
  return (int(var * 100.0) / 100.0);
}

PHTpcSeedFinder::PHTpcSeedFinder()
  : mMaxDistance1(3.0)
  , mTripletAngle1(M_PI / 8)
  , mMinHits1(10)
  , mMaxDistance2(6.0)
  , mTripletAngle2(M_PI / 8)
  , mMinHits2(5)
  , mNThreads(1)
  , mRemoveLoopers(false)
  , mMinLooperRadius(10.0)
  , mMaxLooperRadius(70.0)
{
}

std::vector<kdfinder::TrackCandidate<double>*> PHTpcSeedFinder::findSeeds(TrkrClusterContainer* cluster_map, double B)
{
  LOG_WARN_IF("tracking.PHTpcSeedFnder.findSeeds", !cluster_map) << __FILE__ << "," << __LINE__ << " Can't find node TRKR_CLUSTER";
  std::vector<kdfinder::TrackCandidate<double>*> result;
  if (!cluster_map)
  {
    return result;
  }

  //----- convert clusters to kdhits -----
  std::vector<std::vector<double> > kdhits(PHTpcTrackerUtil::convert_clusters_to_hits(cluster_map));
  LOG_DEBUG("tracking.PHTpcSeedFinder.findSeeds") << "imported " << kdhits.size() << " clusters";

  //----- find unfitted kdtracks -----
  std::vector<std::vector<double> > unused_hits;
  std::vector<std::vector<std::vector<double> > > kdtracks;

  kdtracks = kdfinder::find_tracks_iterative<double>(kdhits, unused_hits,
                                                     mMaxDistance1 /* max distance in cm*/, mTripletAngle1 /* triplet angle */, mMinHits1 /* min hits to keep track */,  // first iteration
                                                     mMaxDistance2, mTripletAngle2, mMinHits2,                                                                           // second iteration params
                                                     mNThreads /* nthreads */,
                                                     false);  // print stats

  LOG_DEBUG("tracking.PHTpcSeedFinder.findSeeds") << "kdtracks: " << kdtracks.size();

  //----- create fitted track candidates -----
  result = kdfinder::get_track_candidates<double>(kdtracks, B);
  LOG_DEBUG("tracking.PHTpcSeedFinder.findSeeds") << "initial kdtrack candidates: " << result.size();

  //----- filter out junk = bad tracks, <25 MeV Pt tracks etc.. -----
  for (auto it = result.begin(); it != result.end();)
  {
    if (!(*it)->isFitted() || (*it)->Pt() < 0.025 || (*it)->Pt() > 200.0)
    {
      (*it)->deleteHits();
      it = result.erase(it);
    }
    else
    {
      ++it;
    }
  }
  LOG_DEBUG("tracking.PHTpcSeedFinder.findSeeds") << "kdtrack candidates after filtering: " << result.size();

  //----- filter out loopers, if enabled -----
  if (mRemoveLoopers)
  {
    for (auto it = result.begin(); it != result.end();)
    {
      kdfinder::Circle<double>* c = (*it)->getCircleFit();
      double x = c->a, y = c->b, r = c->r, xyr = std::sqrt(x * x + y * y);
      if (((xyr - r) > mMinLooperRadius) && ((xyr + r) < mMaxLooperRadius))
      {
        (*it)->deleteHits();
        it = result.erase(it);
      }
      else
      {
        ++it;
      }
    }
  }

  /*
	//----- primitive yet very fast vertex seed finding -----
	// FIXME: use it instead of Rave?
	std::vector< std::pair< std::vector<double>, std::vector<size_t> > > vertices =
    	kdfinder::find_vertex_seeds<double>( result, 0, 0, 0.5, 2.0, 3, false );
    LOG_DEBUG("tracking.PHTpcSeedFinder.findSeeds") << "track candidate vertex seeds: " << vertices.size();
	for( auto vertex: vertices ) {
		std::vector<double> pos = vertex.first;
		std::vector<size_t> ids = vertex.second;
		LOG_DEBUG("tracking.PHTpcSeedFinder.findSeeds") << "vertex tracks: " << ids.size() << ", pos: " << pos[0] << ", " << pos[1] << ", " << pos[2];
	}
	*/

  return result;
}
