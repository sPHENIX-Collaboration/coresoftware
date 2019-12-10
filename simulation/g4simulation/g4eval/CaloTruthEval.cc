#include "CaloTruthEval.h"

#include <phool/getClass.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <utility>

using namespace std;

CaloTruthEval::CaloTruthEval(PHCompositeNode* topNode, const std::string& caloname)
  : _basetrutheval(topNode)
  , _caloname(caloname)
  , _caloid(PHG4HitDefs::get_volume_id(caloname))
  , _truthinfo(nullptr)
  , _g4hits(nullptr)
  , _g4hit_container_id(-1)
  , _strict(false)
  , _verbosity(1)
  , _errors(0)
  , _do_cache(true)
  , _cache_get_shower_energy_deposit()
  , _cache_all_truth_hits_g4shower()
  , _cache_all_truth_hits_g4particle()
  , _cache_get_primary_particle_g4hit()
  , _cache_get_shower_hits_from_primary()
{
  get_node_pointers(topNode);
}

CaloTruthEval::~CaloTruthEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      cout << "CaloTruthEval::~CaloTruthEval() - Error Count: " << _errors << endl;
    }
  }
}

void CaloTruthEval::next_event(PHCompositeNode* topNode)
{
  _cache_get_shower_energy_deposit.clear();
  _cache_all_truth_hits_g4shower.clear();
  _cache_all_truth_hits_g4particle.clear();
  _cache_get_primary_particle_g4hit.clear();
  _cache_get_shower_hits_from_primary.clear();

  _basetrutheval.next_event(topNode);

  get_node_pointers(topNode);
}

bool CaloTruthEval::has_reduced_node_pointers()
{
  if (!_basetrutheval.has_reduced_node_pointers()) return false;

  if (_strict)
    assert(_truthinfo);
  else if (!_truthinfo)
    return false;

  return true;
}

PHG4Shower* CaloTruthEval::get_primary_shower(PHG4Shower* shower)
{
  return _basetrutheval.get_primary_shower(shower);
}

PHG4Shower* CaloTruthEval::get_primary_shower(PHG4Particle* particle)
{
  return _basetrutheval.get_primary_shower(particle);
}

PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Shower* shower)
{
  return _basetrutheval.get_primary_particle(shower);
}

PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Particle* particle)
{
  return _basetrutheval.get_primary_particle(particle);
}

int CaloTruthEval::get_embed(PHG4Particle* particle)
{
  return _basetrutheval.get_embed(particle);
}

PHG4VtxPoint* CaloTruthEval::get_vertex(PHG4Particle* particle)
{
  return _basetrutheval.get_vertex(particle);
}

bool CaloTruthEval::are_same_shower(PHG4Shower* s1, PHG4Shower* s2)
{
  return _basetrutheval.are_same_shower(s1, s2);
}

bool CaloTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2)
{
  return _basetrutheval.are_same_particle(p1, p2);
}

bool CaloTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2)
{
  return _basetrutheval.are_same_vertex(vtx1, vtx2);
}

bool CaloTruthEval::is_primary(PHG4Shower* shower)
{
  return _basetrutheval.is_primary(shower);
}

bool CaloTruthEval::is_primary(PHG4Particle* particle)
{
  return _basetrutheval.is_primary(particle);
}

float CaloTruthEval::get_shower_energy_deposit(PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return NAN;
  }

  if (!is_primary(primary)) return NAN;

  primary = get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, float>::iterator iter =
        _cache_get_shower_energy_deposit.find(primary);
    if (iter != _cache_get_shower_energy_deposit.end())
    {
      return iter->second;
    }
  }

  float shower_e = 0.0;
  PHG4Shower* shower = get_primary_shower(primary);
  if (shower) shower_e = shower->get_edep(_g4hit_container_id);

  if (_do_cache) _cache_get_shower_energy_deposit.insert(make_pair(primary, shower_e));

  return shower_e;
}

float CaloTruthEval::get_shower_eh_ratio(PHG4Particle* primary)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return NAN;
  }

  if (!is_primary(primary)) return NAN;

  PHG4Shower* shower = get_primary_shower(primary);
  if (!shower) return 0.0;

  float ratio = shower->get_eh_ratio(get_caloid());

  return ratio;
}

bool CaloTruthEval::has_full_node_pointers()
{
  if (!_basetrutheval.has_full_node_pointers()) return false;

  if (_strict)
    assert(_truthinfo);
  else if (!_truthinfo)
    return false;

  if (_strict) assert(_g4hits);
  if (!_g4hits) return false;

  return true;
}

std::set<PHG4Hit*> CaloTruthEval::all_truth_hits(PHG4Shower* shower)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }
  else if (!_g4hits)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Shower*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits_g4shower.find(shower);
    if (iter != _cache_all_truth_hits_g4shower.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all g4hits on the shower
  PHG4Shower::HitIdIter iter = shower->find_g4hit_id(_g4hit_container_id);
  if (iter != shower->end_g4hit_id()) return truth_hits;

  for (std::set<PHG4HitDefs::keytype>::iterator jter = iter->second.begin();
       jter != iter->second.end();
       ++jter)
  {
    PHG4Hit* g4hit = _g4hits->findHit(*jter);

    if (_strict)
      assert(g4hit);
    else if (!g4hit)
    {
      ++_errors;
    }

    if (g4hit) truth_hits.insert(g4hit);
  }

  if (_do_cache) _cache_all_truth_hits_g4shower.insert(make_pair(shower, truth_hits));

  return truth_hits;
}

std::set<PHG4Hit*> CaloTruthEval::all_truth_hits(PHG4Particle* particle)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }
  if (!_g4hits)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits_g4particle.find(particle);
    if (iter != _cache_all_truth_hits_g4particle.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits
  for (PHG4HitContainer::ConstIterator g4iter = _g4hits->getHits().first;
       g4iter != _g4hits->getHits().second;
       ++g4iter)
  {
    PHG4Hit* g4hit = g4iter->second;
    if (is_g4hit_from_particle(g4hit, particle)) continue;
    truth_hits.insert(g4hit);
  }

  if (_do_cache) _cache_all_truth_hits_g4particle.insert(make_pair(particle, truth_hits));

  return truth_hits;
}

PHG4Particle* CaloTruthEval::get_parent_particle(PHG4Hit* g4hit)
{
  return _basetrutheval.get_particle(g4hit);
}

PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Hit* g4hit)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, PHG4Particle*>::iterator iter =
        _cache_get_primary_particle_g4hit.find(g4hit);
    if (iter != _cache_get_primary_particle_g4hit.end())
    {
      return iter->second;
    }
  }

  PHG4Particle* primary = _basetrutheval.get_primary_particle(g4hit);

  if (_do_cache) _cache_get_primary_particle_g4hit.insert(make_pair(g4hit, primary));

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
  }

  return primary;
}

bool CaloTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle)
{
  return _basetrutheval.is_g4hit_from_particle(g4hit, particle);
}

std::set<PHG4Hit*> CaloTruthEval::get_shower_hits_from_primary(PHG4Particle* primary)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (!is_primary(primary)) return std::set<PHG4Hit*>();

  primary = get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<PHG4Hit*> >::iterator iter =
        _cache_get_shower_hits_from_primary.find(primary);
    if (iter != _cache_get_shower_hits_from_primary.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  PHG4Shower* shower = get_primary_shower(primary);
  if (shower) truth_hits = all_truth_hits(shower);

  if (_do_cache) _cache_get_shower_hits_from_primary.insert(make_pair(primary, truth_hits));

  return truth_hits;
}

void CaloTruthEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  std::string name = "G4HIT_" + _caloname;
  _g4hits = findNode::getClass<PHG4HitContainer>(topNode, name.c_str());
  _g4hit_container_id = PHG4HitDefs::get_volume_id(name);

  return;
}
