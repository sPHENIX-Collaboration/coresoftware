#include "SvtxTrackEval.h"

#include "SvtxClusterEval.h"

#include <trackbase_historic/SvtxCluster.h>
#include <trackbase_historic/SvtxClusterMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <float.h>
#include <cassert>
#include <set>

using namespace std;

SvtxTrackEval::SvtxTrackEval(PHCompositeNode* topNode)
  : _clustereval(topNode)
  , _trackmap(nullptr)
  , _clustermap(nullptr)
  , _strict(false)
  , _verbosity(1)
  , _errors(0)
  , _do_cache(true)
  , _cache_track_from_cluster_exists(false)
  , _cache_all_truth_hits()
  , _cache_all_truth_particles()
  , _cache_max_truth_particle_by_nclusters()
  , _cache_all_tracks_from_particle()
  , _cache_best_track_from_particle()
  , _cache_all_tracks_from_g4hit()
  , _cache_all_tracks_from_cluster()
  , _cache_best_track_from_cluster()
  , _cache_get_nclusters_contribution()
  , _cache_get_nclusters_contribution_by_layer()
  , _cache_get_nwrongclusters_contribution()
{
  get_node_pointers(topNode);
}

SvtxTrackEval::~SvtxTrackEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      cout << "SvtxTrackEval::~SvtxTrackEval() - Error Count: " << _errors << endl;
    }
  }
}

void SvtxTrackEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_hits.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_nclusters.clear();
  _cache_all_tracks_from_particle.clear();
  _cache_best_track_from_particle.clear();
  _cache_all_tracks_from_g4hit.clear();
  _cache_all_tracks_from_cluster.clear();
  _cache_best_track_from_cluster.clear();
  _cache_get_nclusters_contribution.clear();
  _cache_get_nclusters_contribution_by_layer.clear();
  _cache_get_nwrongclusters_contribution.clear();
  _clustereval.next_event(topNode);

  get_node_pointers(topNode);
}

std::set<PHG4Hit*> SvtxTrackEval::all_truth_hits(SvtxTrack* track)
{
  if (!has_node_pointers()) return std::set<PHG4Hit*>();

  if (_strict)
  {
    assert(track);
  }
  else if (!track)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<SvtxTrack*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(track);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all clusters...
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter)
  {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);

    if (_strict)
    {
      assert(cluster);
    }
    else if (!cluster)
    {
      ++_errors;
      continue;
    }

    std::set<PHG4Hit*> new_hits = _clustereval.all_truth_hits(cluster);

    for (std::set<PHG4Hit*>::iterator jter = new_hits.begin();
         jter != new_hits.end();
         ++jter)
    {
      truth_hits.insert(*jter);
    }
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(track, truth_hits));

  return truth_hits;
}

std::set<PHG4Particle*> SvtxTrackEval::all_truth_particles(SvtxTrack* track)
{
  if (!has_node_pointers()) return std::set<PHG4Particle*>();

  if (_strict)
  {
    assert(track);
  }
  else if (!track)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<SvtxTrack*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(track);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  // loop over all clusters...
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter)
  {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);

    if (_strict)
    {
      assert(cluster);
    }
    else if (!cluster)
    {
      ++_errors;
      continue;
    }

    std::set<PHG4Particle*> new_particles = _clustereval.all_truth_particles(cluster);

    for (std::set<PHG4Particle*>::iterator jter = new_particles.begin();
         jter != new_particles.end();
         ++jter)
    {
      truth_particles.insert(*jter);
    }
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(track, truth_particles));

  return truth_particles;
}

PHG4Particle* SvtxTrackEval::max_truth_particle_by_nclusters(SvtxTrack* track)
{
  if (!has_node_pointers()) return nullptr;

  if (_strict)
  {
    assert(track);
  }
  else if (!track)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<SvtxTrack*, PHG4Particle*>::iterator iter =
        _cache_max_truth_particle_by_nclusters.find(track);
    if (iter != _cache_max_truth_particle_by_nclusters.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> particles = all_truth_particles(track);

  PHG4Particle* max_particle = nullptr;
  unsigned int max_nclusters = 0;

  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
  {
    PHG4Particle* candidate = *iter;
    unsigned int nclusters = get_nclusters_contribution(track, candidate);
    if (nclusters > max_nclusters)
    {
      max_nclusters = nclusters;
      max_particle = candidate;
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_nclusters.insert(make_pair(track, max_particle));

  return max_particle;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Particle* truthparticle)
{
  if (!has_node_pointers()) return std::set<SvtxTrack*>();

  if (_strict)
  {
    assert(truthparticle);
  }
  else if (!truthparticle)
  {
    ++_errors;
    return std::set<SvtxTrack*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<SvtxTrack*> >::iterator iter =
        _cache_all_tracks_from_particle.find(truthparticle);
    if (iter != _cache_all_tracks_from_particle.end())
    {
      return iter->second;
    }
  }

  std::set<SvtxTrack*> tracks;

  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* track = iter->second;

    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
         iter != track->end_clusters();
         ++iter)
    {
      unsigned int cluster_id = *iter;
      SvtxCluster* cluster = _clustermap->get(cluster_id);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      // loop over all particles
      std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster);
      for (std::set<PHG4Particle*>::iterator jter = particles.begin();
           jter != particles.end();
           ++jter)
      {
        PHG4Particle* candidate = *jter;
        if (get_truth_eval()->are_same_particle(candidate, truthparticle))
        {
          tracks.insert(track);
        }
      }
    }
  }

  if (_do_cache) _cache_all_tracks_from_particle.insert(make_pair(truthparticle, tracks));

  return tracks;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Hit* truthhit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<SvtxTrack*>();
  }

  if (_strict)
  {
    assert(truthhit);
  }
  else if (!truthhit)
  {
    ++_errors;
    return std::set<SvtxTrack*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, std::set<SvtxTrack*> >::iterator iter =
        _cache_all_tracks_from_g4hit.find(truthhit);
    if (iter != _cache_all_tracks_from_g4hit.end())
    {
      return iter->second;
    }
  }

  std::set<SvtxTrack*> tracks;

  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* track = iter->second;

    // loop over all clusters
    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
         iter != track->end_clusters();
         ++iter)
    {
      unsigned int cluster_id = *iter;
      SvtxCluster* cluster = _clustermap->get(cluster_id);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      // loop over all hits
      std::set<PHG4Hit*> hits = _clustereval.all_truth_hits(cluster);
      for (std::set<PHG4Hit*>::iterator jter = hits.begin();
           jter != hits.end();
           ++jter)
      {
        PHG4Hit* candidate = *jter;
        // if track id matches argument add to output
        if (candidate->get_trkid() == truthhit->get_trkid())
        {
          tracks.insert(track);
        }
      }
    }
  }

  if (_do_cache) _cache_all_tracks_from_g4hit.insert(make_pair(truthhit, tracks));

  return tracks;
}

SvtxTrack* SvtxTrackEval::best_track_from(PHG4Particle* truthparticle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(truthparticle);
  }
  else if (!truthparticle)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, SvtxTrack*>::iterator iter =
        _cache_best_track_from_particle.find(truthparticle);
    if (iter != _cache_best_track_from_particle.end())
    {
      return iter->second;
    }
  }

  SvtxTrack* best_track = nullptr;
  unsigned int best_count = 0;
  std::set<SvtxTrack*> tracks = all_tracks_from(truthparticle);
  for (std::set<SvtxTrack*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter)
  {
    SvtxTrack* track = *iter;
    unsigned int count = get_nclusters_contribution(track, truthparticle);
    if (count > best_count)
    {
      best_track = track;
      best_count = count;
    }
  }

  if (_do_cache) _cache_best_track_from_particle.insert(make_pair(truthparticle, best_track));

  return best_track;
}

void SvtxTrackEval::create_cache_track_from_cluster()
{
  if (!has_node_pointers())
  {
    ++_errors;
    return;
  }

  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* track = iter->second;

    // loop over all clusters
    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
         iter != track->end_clusters();
         ++iter)
    {
      unsigned int candidate_id = *iter;
      SvtxCluster* candidate = _clustermap->get(candidate_id);

      if (_strict)
      {
        assert(candidate);
      }
      else if (!candidate)
      {
        ++_errors;
        continue;
      }

      //check if cluster has an entry in cache
      std::map<SvtxCluster*, std::set<SvtxTrack*> >::iterator cliter =
          _cache_all_tracks_from_cluster.find(candidate);
      if (cliter != _cache_all_tracks_from_cluster.end())
      {                                //got entry
        cliter->second.insert(track);  //add track to list;
      }
      else
      {
        std::set<SvtxTrack*> tracks;
        tracks.insert(track);
        _cache_all_tracks_from_cluster.insert(make_pair(candidate, tracks));
      }
    }
  }
  _cache_track_from_cluster_exists = true;

  return;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(SvtxCluster* cluster)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<SvtxTrack*>();
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return std::set<SvtxTrack*>();
  }

  std::set<SvtxTrack*> tracks;

  if (_do_cache)
  {
    if (_cache_track_from_cluster_exists == false) create_cache_track_from_cluster();
    std::map<SvtxCluster*, std::set<SvtxTrack*> >::iterator iter =
        _cache_all_tracks_from_cluster.find(cluster);
    if (iter != _cache_all_tracks_from_cluster.end())
    {
      return iter->second;
    }
    else
    {
      return tracks;
    }
  }

  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter)
  {
    SvtxTrack* track = iter->second;

    // loop over all clusters
    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
         iter != track->end_clusters();
         ++iter)
    {
      unsigned int candidate_id = *iter;
      SvtxCluster* candidate = _clustermap->get(candidate_id);

      if (_strict)
      {
        assert(candidate);
      }
      else if (!candidate)
      {
        ++_errors;
        continue;
      }

      if (cluster->get_id() == candidate->get_id())
      {
        tracks.insert(track);
      }
    }
  }

  if (_do_cache) _cache_all_tracks_from_cluster.insert(make_pair(cluster, tracks));

  return tracks;
}

SvtxTrack* SvtxTrackEval::best_track_from(SvtxCluster* cluster)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<SvtxCluster*, SvtxTrack*>::iterator iter =
        _cache_best_track_from_cluster.find(cluster);
    if (iter != _cache_best_track_from_cluster.end())
    {
      return iter->second;
    }
  }

  SvtxTrack* best_track = nullptr;
  float best_quality = FLT_MAX;

  std::set<SvtxTrack*> tracks = all_tracks_from(cluster);
  // loop over all SvtxTracks
  for (std::set<SvtxTrack*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter)
  {
    SvtxTrack* candidate = *iter;

    if (candidate->get_quality() < best_quality)
    {
      best_quality = candidate->get_quality();
      best_track = candidate;
    }
  }

  if (_do_cache) _cache_best_track_from_cluster.insert(make_pair(cluster, best_track));

  return best_track;
}

// overlap calculations
unsigned int SvtxTrackEval::get_nclusters_contribution(SvtxTrack* track, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return 0;
  }

  if (_strict)
  {
    assert(track);
    assert(particle);
  }
  else if (!track || !particle)
  {
    ++_errors;
    return 0;
  }

  calc_cluster_contribution(track, particle);

  std::map<std::pair<SvtxTrack*, PHG4Particle*>, unsigned int>::iterator iter =
      _cache_get_nclusters_contribution.find(make_pair(track, particle));
  if (iter != _cache_get_nclusters_contribution.end())
  {
    return iter->second;
  }

  return 0;
}
unsigned int SvtxTrackEval::get_nwrongclusters_contribution(SvtxTrack* track, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return 0;
  }

  if (_strict)
  {
    assert(track);
    assert(particle);
  }
  else if (!track || !particle)
  {
    ++_errors;
    return 0;
  }

  calc_cluster_contribution(track, particle);

  std::map<std::pair<SvtxTrack*, PHG4Particle*>, unsigned int>::iterator iter =
      _cache_get_nwrongclusters_contribution.find(make_pair(track, particle));
  if (iter != _cache_get_nwrongclusters_contribution.end())
  {
    return iter->second;
  }

  return 0;
}

// overlap calculations
void SvtxTrackEval::calc_cluster_contribution(SvtxTrack* track, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return;
  }

  if (_strict)
  {
    assert(track);
    assert(particle);
  }
  else if (!track || !particle)
  {
    ++_errors;
    return;
  }

  std::map<std::pair<SvtxTrack*, PHG4Particle*>, unsigned int>::iterator iter =
      _cache_get_nclusters_contribution.find(make_pair(track, particle));
  std::map<std::pair<SvtxTrack*, PHG4Particle*>, unsigned int>::iterator witer =
      _cache_get_nwrongclusters_contribution.find(make_pair(track, particle));

  if (iter != _cache_get_nclusters_contribution.end() &&
      witer != _cache_get_nwrongclusters_contribution.end())
  {
    return;
  }

  unsigned int nclusters = 0;
  unsigned int nwrong = 0;
  // loop over all clusters
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter)
  {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);

    if (_strict)
    {
      assert(cluster);
    }
    else if (!cluster)
    {
      ++_errors;
      continue;
    }
    int matched = 0;
    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
         jter != particles.end();
         ++jter)
    {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate, particle))
      {
        ++nclusters;
        matched = 1;
      }
    }
    if (matched == 0) nwrong++;
  }

  _cache_get_nclusters_contribution.insert(make_pair(make_pair(track, particle), nclusters));
  _cache_get_nwrongclusters_contribution.insert(make_pair(make_pair(track, particle), nwrong));

  return;
}

unsigned int SvtxTrackEval::get_nclusters_contribution_by_layer(SvtxTrack* track, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return 0;
  }

  if (_strict)
  {
    assert(track);
    assert(particle);
  }
  else if (!track || !particle)
  {
    ++_errors;
    return 0;
  }

  if (_do_cache)
  {
    std::map<std::pair<SvtxTrack*, PHG4Particle*>, unsigned int>::iterator iter =
        _cache_get_nclusters_contribution_by_layer.find(make_pair(track, particle));
    if (iter != _cache_get_nclusters_contribution_by_layer.end())
    {
      return iter->second;
    }
  }

  unsigned int nclusters_by_layer = 0;

  // loop over all clusters
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter)
  {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);
    unsigned int cluster_layer = cluster->get_layer();

    if (_strict)
    {
      assert(cluster);
    }
    else if (!cluster)
    {
      ++_errors;
      continue;
    }

    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
         jter != particles.end();
         ++jter)
    {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate, particle))
      {
        nclusters_by_layer |= (0x3FFFFFFF & (0x1 << cluster_layer));
      }
    }
  }

  if (_do_cache) _cache_get_nclusters_contribution_by_layer.insert(make_pair(make_pair(track, particle), nclusters_by_layer));

  return nclusters_by_layer;
}

unsigned int SvtxTrackEval::get_layer_range_contribution(SvtxTrack* track, PHG4Particle* particle, unsigned int start_layer, unsigned int end_layer)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return 0;
  }

  if (_strict)
  {
    assert(track);
    assert(particle);
  }
  else if (!track || !particle)
  {
    ++_errors;
    return 0;
  }

  unsigned int nmatches = 0;
  unsigned int nlayers = end_layer - start_layer;

  int layers[nlayers];
  for (unsigned int i = 0; i < nlayers; i++)
  {
    layers[i] = 0;
  }
  // loop over all clusters
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter)
  {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);
    unsigned int cluster_layer = cluster->get_layer();
    if (cluster_layer >= end_layer) continue;
    if (cluster_layer < start_layer) continue;

    if (_strict)
    {
      assert(cluster);
    }
    else if (!cluster)
    {
      ++_errors;
      continue;
    }

    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
         jter != particles.end();
         ++jter)
    {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate, particle))
      {
        //	nmatches |= (0x3FFFFFFF & (0x1 << cluster_layer));
        layers[cluster_layer - start_layer] = 1;
      }
    }
  }
  for (unsigned int i = 0; i < nlayers; i++)
  {
    if (layers[i] == 1) nmatches++;
  }

  return nmatches;
}

void SvtxTrackEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  _clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");

  return;
}

bool SvtxTrackEval::has_node_pointers()
{
  // need things off of the DST...
  if (_strict)
    assert(_trackmap);
  else if (!_trackmap)
    return false;

  if (_strict)
    assert(_clustermap);
  else if (!_clustermap)
    return false;

  return true;
}
