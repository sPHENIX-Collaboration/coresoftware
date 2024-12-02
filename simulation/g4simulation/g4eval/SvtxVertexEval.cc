#include "SvtxVertexEval.h"

#include "SvtxTrackEval.h"
#include "SvtxTruthEval.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <globalvertex/Vertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>

#include <cassert>
#include <iostream>
#include <set>

class SvtxTrack;

SvtxVertexEval::SvtxVertexEval(PHCompositeNode* topNode)
  : _trackeval(topNode)
{
  set_track_nodename("SvtxTrackMap");
}

SvtxVertexEval::~SvtxVertexEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      std::cout << "SvtxVertexEval::~SvtxVertexEval() - Error Count: " << _errors << std::endl;
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

std::set<PHG4Particle*> SvtxVertexEval::all_truth_particles(const Vertex* vertex)
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
    std::map<const Vertex*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(vertex);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> all_particles;

  // loop over all tracks on vertex
  for (Vertex::TrackIter iter = vertex->begin_tracks();
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

  if (_do_cache)
  {
    _cache_all_truth_particles.insert(std::make_pair(vertex, all_particles));
  }

  return all_particles;
}

std::set<PHG4VtxPoint*> SvtxVertexEval::all_truth_points(const Vertex* vertex)
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
    std::map<const Vertex*, std::set<PHG4VtxPoint*> >::iterator iter =
        _cache_all_truth_points.find(vertex);
    if (iter != _cache_all_truth_points.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4VtxPoint*> points;

  std::set<PHG4Particle*> particles = all_truth_particles(vertex);
  for (auto particle : particles)
  {
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

  if (_do_cache)
  {
    _cache_all_truth_points.insert(std::make_pair(vertex, points));
  }

  return points;
}

PHG4VtxPoint* SvtxVertexEval::max_truth_point_by_ntracks(const Vertex* vertex)
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
    std::map<const Vertex*, PHG4VtxPoint*>::iterator iter =
        _cache_max_truth_point_by_ntracks.find(vertex);
    if (iter != _cache_max_truth_point_by_ntracks.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4VtxPoint*> points = all_truth_points(vertex);

  PHG4VtxPoint* max_point = nullptr;
  unsigned int max_ntracks = 0;

  for (auto candidate : points)
  {
    unsigned int ntracks = get_ntracks_contribution(vertex, candidate);
    if (ntracks > max_ntracks)
    {
      max_ntracks = ntracks;
      max_point = candidate;
    }
  }

  if (_do_cache)
  {
    _cache_max_truth_point_by_ntracks.insert(std::make_pair(vertex, max_point));
  }

  return max_point;
}

std::set<const Vertex*> SvtxVertexEval::all_vertexes_from(PHG4VtxPoint* truthpoint)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<const Vertex*>();
  }

  if (_strict)
  {
    assert(truthpoint);
  }
  else if (!truthpoint)
  {
    ++_errors;
    return std::set<const Vertex*>();
  }

  if (_do_cache)
  {
    std::map<PHG4VtxPoint*, std::set<const Vertex*> >::iterator iter =
        _cache_all_vertexes_from_point.find(truthpoint);
    if (iter != _cache_all_vertexes_from_point.end())
    {
      return iter->second;
    }
  }

  std::set<const Vertex*> all_vertexes;

  // loop over all vertexes on node

  for (auto& iter : *_vertexmap)
  {
    const Vertex* vertex = iter.second;
    std::set<PHG4VtxPoint*> points = all_truth_points(vertex);
    for (auto point : points)
    {
      if (get_truth_eval()->are_same_vertex(point, truthpoint))
      {
        all_vertexes.insert(vertex);
      }
    }
  }

  if (_do_cache)
  {
    _cache_all_vertexes_from_point.insert(std::make_pair(truthpoint, all_vertexes));
  }

  return all_vertexes;
}

const Vertex* SvtxVertexEval::best_vertex_from(PHG4VtxPoint* truthpoint)
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
    std::map<PHG4VtxPoint*, const Vertex*>::iterator iter =
        _cache_best_vertex_from_point.find(truthpoint);
    if (iter != _cache_best_vertex_from_point.end())
    {
      return iter->second;
    }
  }

 
  unsigned int best_count = 0;
  std::set<const Vertex*> tracks = all_vertexes_from(truthpoint);
  
  std::set<const Vertex*>::iterator best_vertex_iter = tracks.begin();
  for (auto it = tracks.begin(); it != tracks.end(); ++it)
  {
    const Vertex* vertex = *it;
    unsigned int count = get_ntracks_contribution(vertex, truthpoint);
    if (count > best_count)
    {
      best_vertex_iter = it;
      best_count = count;
    }
  }

  const Vertex* the_best = *best_vertex_iter;

  if (_do_cache)
  {
    _cache_best_vertex_from_point.insert(std::make_pair(truthpoint, the_best));
  }

  return the_best;
}

// overlap calculations
unsigned int SvtxVertexEval::get_ntracks_contribution(const Vertex* vertex, PHG4VtxPoint* truthpoint)
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
    std::map<std::pair<const Vertex*, PHG4VtxPoint*>, unsigned int>::iterator iter =
        _cache_get_ntracks_contribution.find(std::make_pair(vertex, truthpoint));
    if (iter != _cache_get_ntracks_contribution.end())
    {
      return iter->second;
    }
  }

  unsigned int ntracks = 0;

  for (Vertex::TrackIter iter = vertex->begin_tracks();
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

  if (_do_cache)
  {
    _cache_get_ntracks_contribution.insert(std::make_pair(std::make_pair(vertex, truthpoint), ntracks));
  }

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
  {
    std::cout << PHWHERE << "Did not find_vertexmap on node tree" << std::endl;
  }

  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_TrackNodeName);

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return;
}

bool SvtxVertexEval::has_node_pointers()
{
  if (_strict)
  {
    assert(_vertexmap);
  }
  else if (!_vertexmap)
  {
    std::cout << PHWHERE << " did not find _vertexmap " << std::endl;
    return false;
  }

  if (_strict)
  {
    assert(_trackmap);
  }
  else if (!_trackmap)
  {
    std::cout << PHWHERE << " did not find _trackmap " << std::endl;
    return false;
  }

  if (_strict)
  {
    assert(_truthinfo);
  }
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
