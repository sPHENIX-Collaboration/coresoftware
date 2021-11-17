#include "SvtxVertexEval.h"

#include "SvtxTrackEval.h"
#include "SvtxTruthEval.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>

#include <cassert>
#include <iostream>
#include <set>

class SvtxTrack;

using namespace std;

SvtxVertexEval::SvtxVertexEval(PHCompositeNode* topNode)
  : _trackeval(topNode)
  , _vertexmap(nullptr)
  , _trackmap(nullptr)
  , _truthinfo(nullptr)
  , _strict(false)
  , _verbosity(0)
  , _errors(0)
  , _do_cache(true)
{
  set_track_nodename("SvtxTrackMap");
}

SvtxVertexEval::~SvtxVertexEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      cout << "SvtxVertexEval::~SvtxVertexEval() - Error Count: " << _errors << endl;
    }
  }
}

void SvtxVertexEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_particles.clear();
  _cache_all_truth_points.clear();
  _cache_max_truth_point_by_ntracks.clear();
  _cache_all_vertexes_from_point.clear();
  _cache_best_vertex_from_point.clear();
  _cache_get_ntracks_contribution.clear();

  _trackeval.next_event(topNode);

  get_node_pointers(topNode);
}

std::set<PHG4Particle*> SvtxVertexEval::all_truth_particles(SvtxVertex* vertex)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(vertex);
  }
  else if (!vertex)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<SvtxVertex*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(vertex);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> all_particles;

  // loop over all tracks on vertex
  for (SvtxVertex::TrackIter iter = vertex->begin_tracks();
       iter != vertex->end_tracks();
       ++iter)
  {
    SvtxTrack* track = _trackmap->get(*iter);

    if (_strict)
    {
      assert(track);
    }
    else if (!track)
    {
      ++_errors;
      continue;
    }

    PHG4Particle* max_particle = _trackeval.max_truth_particle_by_nclusters(track);

    if (_strict)
    {
      assert(max_particle);
    }
    else if (!max_particle)
    {
      ++_errors;
      continue;
    }

    all_particles.insert(max_particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(vertex, all_particles));

  return all_particles;
}

std::set<PHG4VtxPoint*> SvtxVertexEval::all_truth_points(SvtxVertex* vertex)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<PHG4VtxPoint*>();
  }

  if (_strict)
  {
    assert(vertex);
  }
  else if (!vertex)
  {
    ++_errors;
    return std::set<PHG4VtxPoint*>();
  }

  if (_do_cache)
  {
    std::map<SvtxVertex*, std::set<PHG4VtxPoint*> >::iterator iter =
        _cache_all_truth_points.find(vertex);
    if (iter != _cache_all_truth_points.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4VtxPoint*> points;

  std::set<PHG4Particle*> particles = all_truth_particles(vertex);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;
    PHG4VtxPoint* point = get_truth_eval()->get_vertex(particle);

    if (_strict)
    {
      assert(point);
    }
    else if (!point)
    {
      ++_errors;
      continue;
    }

    points.insert(point);
  }

  if (_do_cache) _cache_all_truth_points.insert(make_pair(vertex, points));

  return points;
}

PHG4VtxPoint* SvtxVertexEval::max_truth_point_by_ntracks(SvtxVertex* vertex)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(vertex);
  }
  else if (!vertex)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<SvtxVertex*, PHG4VtxPoint*>::iterator iter =
        _cache_max_truth_point_by_ntracks.find(vertex);
    if (iter != _cache_max_truth_point_by_ntracks.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4VtxPoint*> points = all_truth_points(vertex);

  PHG4VtxPoint* max_point = nullptr;
  unsigned int max_ntracks = 0;

  for (std::set<PHG4VtxPoint*>::iterator iter = points.begin();
       iter != points.end();
       ++iter)
  {
    PHG4VtxPoint* candidate = *iter;
    unsigned int ntracks = get_ntracks_contribution(vertex, candidate);
    if (ntracks > max_ntracks)
    {
      max_ntracks = ntracks;
      max_point = candidate;
    }
  }

  if (_do_cache) _cache_max_truth_point_by_ntracks.insert(make_pair(vertex, max_point));

  return max_point;
}

std::set<SvtxVertex*> SvtxVertexEval::all_vertexes_from(PHG4VtxPoint* truthpoint)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<SvtxVertex*>();
  }

  if (_strict)
  {
    assert(truthpoint);
  }
  else if (!truthpoint)
  {
    ++_errors;
    return std::set<SvtxVertex*>();
  }

  if (_do_cache)
  {
    std::map<PHG4VtxPoint*, std::set<SvtxVertex*> >::iterator iter =
        _cache_all_vertexes_from_point.find(truthpoint);
    if (iter != _cache_all_vertexes_from_point.end())
    {
      return iter->second;
    }
  }

  std::set<SvtxVertex*> all_vertexes;

  // loop over all vertexes on node

  for (SvtxVertexMap::Iter iter = _vertexmap->begin();
       iter != _vertexmap->end();
       ++iter)
  {
    SvtxVertex* vertex = iter->second;
    std::set<PHG4VtxPoint*> points = all_truth_points(vertex);
    for (std::set<PHG4VtxPoint*>::iterator jter = points.begin();
         jter != points.end();
         ++jter)
    {
      PHG4VtxPoint* point = *jter;
      if (get_truth_eval()->are_same_vertex(point, truthpoint))
      {
        all_vertexes.insert(vertex);
      }
    }
  }

  if (_do_cache) _cache_all_vertexes_from_point.insert(make_pair(truthpoint, all_vertexes));

  return all_vertexes;
}

SvtxVertex* SvtxVertexEval::best_vertex_from(PHG4VtxPoint* truthpoint)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(truthpoint);
  }
  else if (!truthpoint)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4VtxPoint*, SvtxVertex*>::iterator iter =
        _cache_best_vertex_from_point.find(truthpoint);
    if (iter != _cache_best_vertex_from_point.end())
    {
      return iter->second;
    }
  }

  SvtxVertex* best_vertex = nullptr;
  unsigned int best_count = 0;
  std::set<SvtxVertex*> tracks = all_vertexes_from(truthpoint);
  for (std::set<SvtxVertex*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter)
  {
    SvtxVertex* vertex = *iter;
    unsigned int count = get_ntracks_contribution(vertex, truthpoint);
    if (count > best_count)
    {
      best_vertex = vertex;
      best_count = count;
    }
  }

  if (_do_cache) _cache_best_vertex_from_point.insert(make_pair(truthpoint, best_vertex));

  return best_vertex;
}

// overlap calculations
unsigned int SvtxVertexEval::get_ntracks_contribution(SvtxVertex* vertex, PHG4VtxPoint* truthpoint)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return 0;
  }

  if (_strict)
  {
    assert(vertex);
    assert(truthpoint);
  }
  else if (!vertex || !truthpoint)
  {
    ++_errors;
    return 0;
  }

  if (_do_cache)
  {
    std::map<std::pair<SvtxVertex*, PHG4VtxPoint*>, unsigned int>::iterator iter =
        _cache_get_ntracks_contribution.find(make_pair(vertex, truthpoint));
    if (iter != _cache_get_ntracks_contribution.end())
    {
      return iter->second;
    }
  }

  unsigned int ntracks = 0;

  for (SvtxVertex::TrackIter iter = vertex->begin_tracks();
       iter != vertex->end_tracks();
       ++iter)
  {
    SvtxTrack* track = _trackmap->get(*iter);
    PHG4Particle* particle = _trackeval.max_truth_particle_by_nclusters(track);

    PHG4VtxPoint* candidate = get_truth_eval()->get_vertex(particle);

    if (_strict)
    {
      assert(candidate);
    }
    else if (!candidate)
    {
      ++_errors;
      continue;
    }

    if (get_truth_eval()->are_same_vertex(candidate, truthpoint))
    {
      ++ntracks;
    }
  }

  if (_do_cache) _cache_get_ntracks_contribution.insert(make_pair(make_pair(vertex, truthpoint), ntracks));

  return ntracks;
}

void SvtxVertexEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off the DST...

  if (_use_initial_vertex)
  {
    _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");  // always there, initial vertices
  }
  else if (_use_genfit_vertex)
  {
    _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");  // Rave vertices
  }
  else
  {
    _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapActs");  // Acts vertices
  }
  if (!_vertexmap)
    std::cout << PHWHERE << "Did not find_vertexmap on node tree" << endl;

  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_TrackNodeName);

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return;
}

bool SvtxVertexEval::has_node_pointers()
{
  if (_strict)
    assert(_vertexmap);
  else if (!_vertexmap)
  {
    std::cout << PHWHERE << " did not find _vertexmap " << std::endl;
    return false;
  }

  if (_strict)
    assert(_trackmap);
  else if (!_trackmap)
  {
    std::cout << PHWHERE << " did not find _trackmap " << std::endl;
    return false;
  }

  if (_strict)
    assert(_truthinfo);
  else if (!_truthinfo)
  {
    std::cout << PHWHERE << " did not find _truthinfo " << std::endl;
    return false;
  }

  return true;
}

void SvtxVertexEval::set_track_nodename(const std::string& name)
{
  m_TrackNodeName = name;
  _trackeval.set_track_nodename(name);
}
