#include "SvtxTrackEval.h"

#include "SvtxClusterEval.h"
#include "SvtxTruthEval.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <trackbase/TrkrDefs.h>  // for cluskey, getLayer

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxPHG4ParticleMap.h>
#include <trackbase_historic/PHG4ParticleSvtxMap.h>


#include <phool/getClass.h>

#include <cassert>
#include <cfloat>
#include <iostream>
#include <set>

using namespace std;

SvtxTrackEval::SvtxTrackEval(PHCompositeNode* topNode)
  : _clustereval(topNode)
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
  std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
 
  // loop over all clusters...
  for (const auto& cluster_key : cluster_keys)
    {
    //    if (_strict)
    //    {
    //      assert(cluster_key);
    //    }
    //    else if (!cluster_key)
    //    {
    //      ++_errors;
    //      continue;
    //    }

    std::set<PHG4Hit*> new_hits = _clustereval.all_truth_hits(cluster_key);

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

  if(_recoTruthMap && _recoTruthMap->processed()) 
    {
      SvtxPHG4ParticleMap::WeightedTruthTrackMap map = _recoTruthMap->get(track->get_id());
      std::set<PHG4Particle*> returnset;
      
      for(const auto& [weight, truthTrackSet] : map)
	{
	  for(const int& g4partid : truthTrackSet)
	    {
	      returnset.insert(_truthinfo->GetParticle(g4partid));
	    }
	}
      return returnset;
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
  SvtxTrack_FastSim * fastsim_track = dynamic_cast<SvtxTrack_FastSim * >(track);

  if (fastsim_track)
  {
    // exception for fast sim track
    unsigned int track_id = fastsim_track -> get_truth_track_id();
    truth_particles.insert(get_truth_eval()->get_particle(track_id));
  }
  else{                
    // loop over all clusters...
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
    for (const auto& cluster_key : cluster_keys)
    {
      //    if (_strict)
      //    {
      //      assert(cluster_key);
      //    }
      //    else if (!cluster_key)
      //    {
      //      ++_errors;
      //      continue;
      //    }

      std::set<PHG4Particle*> new_particles = _clustereval.all_truth_particles(cluster_key);

      for (std::set<PHG4Particle*>::iterator jter = new_particles.begin();
          jter != new_particles.end();
          ++jter)
      {
        truth_particles.insert(*jter);
      }
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

  if(_recoTruthMap && _recoTruthMap->processed())
    {
      const SvtxPHG4ParticleMap::WeightedTruthTrackMap map = _recoTruthMap->get(track->get_id());
      if (map.size() == 0) return nullptr;
      auto itr = map.end();
      --itr;
      std::set<int> bestPartSet = itr->second;
      int bestpart = *bestPartSet.begin();
      return _truthinfo->GetParticle(bestpart);
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
  
  SvtxTrack_FastSim * fastsim_track = dynamic_cast<SvtxTrack_FastSim * >(track);
  if (fastsim_track)
  {
    // exception for fast sim track
    unsigned int track_id = fastsim_track -> get_truth_track_id();
    max_particle = get_truth_eval()->get_particle(track_id);
  }
  else
  {
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

  if(_truthRecoMap && _truthRecoMap->processed())
    {
      std::set<SvtxTrack*> returnset;
 
      PHG4ParticleSvtxMap::WeightedRecoTrackMap map = _truthRecoMap->get(truthparticle->get_track_id());    
      
      for(const auto& [weight, recoTrackSet] : map)
	{
	  for(const unsigned int& trackid : recoTrackSet)
	    {
	      returnset.insert(_trackmap->get(trackid));
	    }
	}
      return returnset;
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
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
    for (const auto& cluster_key : cluster_keys)
    {
      // remove this check as cluster key = 0 is MVTX layer 0 cluster #0.
      //      if (_strict)
      //      {
      //        assert(cluster_key);
      //      }
      //      else if (!cluster_key)
      //      {
      //        ++_errors;
      //        continue;
      //      }

      // loop over all particles
      std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster_key);
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
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
    // loop over all clusters
    for (const auto& cluster_key : cluster_keys)
    {

      //      if (_strict)
      //      {
      //        assert(cluster_key);
      //      }
      //      else if (!cluster_key)
      //      {
      //        ++_errors;
      //        continue;
      //      }

      // loop over all hits
      std::set<PHG4Hit*> hits = _clustereval.all_truth_hits(cluster_key);
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
  
  if(_truthRecoMap && _truthRecoMap->processed())
    {
      const PHG4ParticleSvtxMap::WeightedRecoTrackMap map = _truthRecoMap->get(truthparticle->get_track_id());
      /// No reco tracks found
      if(map.size() == 0)
	{ return nullptr; }
      auto itr = map.end();
      --itr;
      std::set<unsigned int> bestPartSet = itr->second;
      int bestpart = *bestPartSet.begin();
      return _trackmap->get(bestpart);
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
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate_key : cluster_keys)
    {
      // unsigned int cluster_layer = TrkrDefs::getLayer(candidate_key);
      //      if (_strict)
      //      {
      //        assert(candidate_key);
      //      }
      //      else if (!candidate_key)
      //      {
      //        ++_errors;
      //        continue;
      //      }

      //check if cluster has an entry in cache
      std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator cliter =
          _cache_all_tracks_from_cluster.find(candidate_key);
      if (cliter != _cache_all_tracks_from_cluster.end())
      {                                //got entry
        cliter->second.insert(track);  //add track to list;
      }
      else
      {
        std::set<SvtxTrack*> tracks;
        tracks.insert(track);
        _cache_all_tracks_from_cluster.insert(make_pair(candidate_key, tracks));
      }
    }
  }
  _cache_track_from_cluster_exists = true;

  return;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<SvtxTrack*>();
  }

  //  if (_strict)
  //  {
  //    assert(cluster_key);
  //  }
  //  else if (!cluster_key)
  //  {
  //    ++_errors;
  //    return std::set<SvtxTrack*>();
  //  }

  std::set<SvtxTrack*> tracks;

  if (_do_cache)
  {
    if (_cache_track_from_cluster_exists == false) create_cache_track_from_cluster();
    std::map<TrkrDefs::cluskey, std::set<SvtxTrack*> >::iterator iter =
        _cache_all_tracks_from_cluster.find(cluster_key);
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
    std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);

    // loop over all clusters
    for (const auto& candidate : cluster_keys)
    {

      //      if (_strict)
      //      {
      //        assert(candidate);
      //      }
      //      else if (!candidate)
      //      {
      //        ++_errors;
      //        continue;
      //      }

      if (cluster_key == candidate)
      {
        tracks.insert(track);
      }
    }
  }

  if (_do_cache) _cache_all_tracks_from_cluster.insert(make_pair(cluster_key, tracks));

  return tracks;
}

SvtxTrack* SvtxTrackEval::best_track_from(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  //  if (_strict)
  //  {
  //    assert(cluster_key);
  //  }
  //  else if (!cluster_key)
  //  {
  //    ++_errors;
  //    return nullptr;
  //  }

  if (_do_cache)
  {
    std::map<TrkrDefs::cluskey, SvtxTrack*>::iterator iter =
        _cache_best_track_from_cluster.find(cluster_key);
    if (iter != _cache_best_track_from_cluster.end())
    {
      return iter->second;
    }
  }

  SvtxTrack* best_track = nullptr;
  float best_quality = FLT_MAX;

  std::set<SvtxTrack*> tracks = all_tracks_from(cluster_key);
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

  if (_do_cache) _cache_best_track_from_cluster.insert(make_pair(cluster_key, best_track));
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
  std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
  for (const auto& cluster_key : cluster_keys)
  {
    //    if (_strict)
    //    {
    //      assert(cluster_key);
    //    }
    //    else if (!cluster_key)
    //    {
    //      ++_errors;
    //      continue;
    //    }
    int matched = 0;
    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster_key);
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
  int layer_occupied[100];
  for (int i = 0; i < 100; i++) layer_occupied[i] = 0;

  // loop over all clusters
  std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
  for (const auto& cluster_key : cluster_keys)
  {
    unsigned int cluster_layer = TrkrDefs::getLayer(cluster_key);

    //    if (_strict)
    //    {
    //      assert(cluster_key);
    //    }
    //    else if (!cluster_key)
    //    {
    //      ++_errors;
    //      continue;
    //    }

    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster_key);

    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
         jter != particles.end();
         ++jter)
    {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate, particle))
      {
        layer_occupied[cluster_layer]++;
      }
    }
  }
  for (int i = 0; i < 100; i++)
  {
    if (layer_occupied[i] > 0) nclusters_by_layer++;
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
  std::vector<TrkrDefs::cluskey> cluster_keys = get_track_ckeys(track);
  for (const auto& cluster_key : cluster_keys)
  {
    unsigned int cluster_layer = TrkrDefs::getLayer(cluster_key);
    if (cluster_layer >= end_layer) continue;
    if (cluster_layer < start_layer) continue;

    //    if (_strict)
    //    {
    //      assert(cluster_key);
    //    }
    //    else if (!cluster_key)
    //    {
    //      ++_errors;
    //      continue;
    //    }

    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster_key);
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
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_TrackNodeName);

  _truthRecoMap = findNode::getClass<PHG4ParticleSvtxMap>(topNode, "PHG4ParticleSvtxMap");

  _recoTruthMap = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
 
  return;
}

bool SvtxTrackEval::has_node_pointers()
{
  // need things off of the DST...
  if (_strict)
    assert(_trackmap);
  else if (!_trackmap)
    return false;
  return true;
}

std::vector<TrkrDefs::cluskey> SvtxTrackEval::get_track_ckeys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> cluster_keys;
  TrackSeed *tpcseed = track->get_tpc_seed();
  TrackSeed *silseed = track->get_silicon_seed();
  if(silseed)
    {
      for(auto iter = silseed->begin_cluster_keys();
	  iter!= silseed->end_cluster_keys();
	  ++iter)
	{ cluster_keys.push_back(*iter); }
    }
  if(tpcseed)
    {
      for(auto iter = tpcseed->begin_cluster_keys();
          iter!= tpcseed->end_cluster_keys();
          ++iter)
        { cluster_keys.push_back(*iter); }
    }
  
  return cluster_keys;
}
