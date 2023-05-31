#include "SvtxHitEval.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>  // for operator<<, endl, basic_...
#include <map>
#include <set>

class TrkrHit;

using namespace std;

SvtxHitEval::SvtxHitEval(PHCompositeNode* topNode)
  : _trutheval(topNode)
  //  , _g4cells_svtx(nullptr)
  // , _g4cells_tracker(nullptr)
  //, _g4cells_maps(nullptr)
  , _g4hits_tpc(nullptr)
  , _g4hits_intt(nullptr)
  , _g4hits_mvtx(nullptr)
  , _g4hits_mms(nullptr)
  , _truthinfo(nullptr)
  , _strict(false)
  , _verbosity(0)
  , _errors(0)
  , _do_cache(true)
  , _cache_all_truth_hits()
  , _cache_max_truth_hit_by_energy()
  , _cache_all_truth_particles()
  , _cache_max_truth_particle_by_energy()
  , _cache_all_hits_from_particle()
  , _cache_all_hits_from_g4hit()
  , _cache_best_hit_from_g4hit()
  , _cache_get_energy_contribution_g4particle()
  , _cache_get_energy_contribution_g4hit()
{
  get_node_pointers(topNode);
}

SvtxHitEval::~SvtxHitEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      cout << "SvtxHitEval::~SvtxHitEval() - Error Count: " << _errors << endl;
    }
  }
}

void SvtxHitEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_hits.clear();
  _cache_max_truth_hit_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_energy.clear();
  _cache_all_hits_from_particle.clear();
  _cache_all_hits_from_g4hit.clear();
  _cache_best_hit_from_g4hit.clear();
  _cache_get_energy_contribution_g4particle.clear();
  _cache_get_energy_contribution_g4hit.clear();

  _trutheval.next_event(topNode);

  get_node_pointers(topNode);
}

std::set<PHG4Hit*> SvtxHitEval::all_truth_hits(const TrkrDefs::hitsetkey hitset_key, TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Hit*>();
  }

  const std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> key(hitset_key, hit_key);

  if (_do_cache)
  {
    const auto &  iter =
        _cache_all_truth_hits.find(key);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  const uint8_t trkrid = TrkrDefs::getTrkrId(hitset_key);

  std::set<PHG4Hit*> truth_hits;

  assert(_hit_truth_map);

  std::multimap<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;
  _hit_truth_map->getG4Hits(hitset_key, hit_key, temp_map);  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
  for (std::multimap<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter = temp_map.begin();
       htiter != temp_map.end(); ++htiter)
  {
    // extract the g4 hit key here and add the g4hit to the set
    PHG4HitDefs::keytype g4hitkey = htiter->second.second;
    // cout << "           hitkey " << hitkey <<  " g4hitkey " << g4hitkey << endl;
    PHG4Hit* g4hit = nullptr;
    switch (trkrid)
    {
    case TrkrDefs::tpcId:
      g4hit = _g4hits_tpc->findHit(g4hitkey);
      break;
    case TrkrDefs::inttId:
      g4hit = _g4hits_intt->findHit(g4hitkey);
      break;
    case TrkrDefs::mvtxId:
      g4hit = _g4hits_mvtx->findHit(g4hitkey);
      break;
    case TrkrDefs::micromegasId:
      g4hit = _g4hits_mms->findHit(g4hitkey);
      break;
    default:
      break;
    }
    // fill output set
    if (g4hit) truth_hits.insert(g4hit);
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(key, truth_hits));

  return truth_hits;
}

PHG4Hit* SvtxHitEval::max_truth_hit_by_energy(const TrkrDefs::hitsetkey hitset_key, TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  const std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> key(hitset_key, hit_key);

  if (_do_cache)
  {
    const auto &  iter =
        _cache_max_truth_hit_by_energy.find(key);
    if (iter != _cache_max_truth_hit_by_energy.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> hits = all_truth_hits(hitset_key, hit_key);
  PHG4Hit* max_hit = nullptr;
  float max_e = FLT_MAX * -1.0;
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter)
  {
    PHG4Hit* hit = *iter;
    if (hit->get_edep() > max_e)
    {
      max_e = hit->get_edep();
      max_hit = hit;
    }
  }

  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(key, max_hit));

  return max_hit;
}

std::set<PHG4Particle*> SvtxHitEval::all_truth_particles(const TrkrDefs::hitsetkey hitset_key, TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Particle*>();
  }

  const std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> key(hitset_key, hit_key);

  if (_do_cache)
  {
    const auto &  iter =
        _cache_all_truth_particles.find(key);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  std::set<PHG4Hit*> g4hits = all_truth_hits(hitset_key, hit_key);

  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
  {
    PHG4Hit* g4hit = *iter;
    PHG4Particle* particle = get_truth_eval()->get_particle(g4hit);

    if (_strict)
      assert(particle);
    else if (!particle)
    {
      ++_errors;
      if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
      continue;
    }

    truth_particles.insert(particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(key, truth_particles));

  return truth_particles;
}

PHG4Particle* SvtxHitEval::max_truth_particle_by_energy(const TrkrDefs::hitsetkey hitset_key, TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  const std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> key(hitset_key, hit_key);
  if (_do_cache)
  {
    const auto &  iter =
        _cache_max_truth_particle_by_energy.find(key);
    if (iter != _cache_max_truth_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  // loop over all particles associated with this hit and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(hitset_key, hit_key);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(hitset_key, hit_key, particle);
    if (e > max_e)
    {
      max_e = e;
      max_particle = particle;
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(key, max_particle));

  return max_particle;
}

//std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> > SvtxHitEval::all_hits_from(PHG4Particle* g4particle)
//{
//  if (!has_node_pointers())
//  {
//    ++_errors;
//    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
//    return std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> >();
//  }
//
//  if (_strict)
//  {
//    assert(g4particle);
//  }
//  else if (!g4particle)
//  {
//    ++_errors;
//    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
//    return std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> >();
//  }
//
//  if (_do_cache)
//  {
//    const auto & iter =
//        _cache_all_hits_from_particle.find(g4particle);
//    if (iter != _cache_all_hits_from_particle.end())
//    {
//      return iter->second;
//    }
//  }
//
//  std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> > hits;
//
//  assert(_hit_truth_map);
//
//  for (const auto& [trkrID, _hitmap] : _hitmaps)
//  {
//    // loop over all the hits
//    TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
//    for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
//         iter != all_hitsets.second;
//         ++iter)
//    {
//      TrkrHitSet::ConstRange range = iter->second->getHits();
//      for (TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
//      {
//        TrkrDefs::hitkey hit_key = hitr->first;
//
//        // loop over all truth hits connected to this hit
//        std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key, trkrID);
//        for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
//             jter != g4hits.end();
//             ++jter)
//        {
//          PHG4Hit* candidate = *jter;
//          PHG4Particle* particle = _truthinfo->GetParticle(candidate->get_trkid());
//          if (g4particle->get_track_id() == particle->get_track_id())
//          {
//            hits.insert(make_pair(iter->first, hit_key));
//          }
//        }
//      }
//    }
//  }
//
//  if (_do_cache) _cache_all_hits_from_particle.insert(make_pair(g4particle, hits));
//
//  return hits;
//}

//std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> > SvtxHitEval::all_hits_from(PHG4Hit* g4hit)
//{
//  if (!has_node_pointers())
//  {
//    ++_errors;
//    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
//    return std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> >();
//  }
//
//  if (_strict)
//  {
//    assert(g4hit);
//  }
//  else if (!g4hit)
//  {
//    ++_errors;
//    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
//    return std::set<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> >();
//  }
//
//  if (_do_cache)
//  {
//    const auto & iter =
//        _cache_all_hits_from_g4hit.find(g4hit);
//    if (iter != _cache_all_hits_from_g4hit.end())
//    {
//      return iter->second;
//    }
//  }
//
//  std::set<TrkrDefs::hitkey> hits;
//
//  unsigned int hit_layer = g4hit->get_layer();
//
//  // loop over all the hits
//  for (const auto& [trkrID, _hitmap] : _hitmaps)
//  {
//    TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
//    for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
//         iter != all_hitsets.second;
//         ++iter)
//    {
//      TrkrHitSet::ConstRange range = iter->second->getHits();
//      for (TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
//      {
//        TrkrDefs::hitkey hit_key = hitr->first;
//
//        if (TrkrDefs::getLayer(hit_key) != hit_layer) continue;
//
//        // loop over all truth hits connected to this hit
//        std::set<PHG4Hit*> g4hits = all_truth_hits(iter->first, hit_key);
//        for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
//             jter != g4hits.end();
//             ++jter)
//        {
//          PHG4Hit* candidate = *jter;
//          if (candidate->get_hit_id() == g4hit->get_hit_id())
//          {
//            hits.insert(make_pair(iter->first, hit_key));
//          }
//        }
//      }
//    }
//  }
//
//  if (_do_cache) _cache_all_hits_from_g4hit.insert(make_pair(g4hit, hits));
//
//  return hits;
//}

// overlap calculations
float SvtxHitEval::get_energy_contribution(const TrkrDefs::hitsetkey hitset_key, TrkrDefs::hitkey hit_key, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }

  if (_strict)
  {
    assert(hit_key);
    assert(particle);
  }
  else if (!hit_key || !particle)
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }

  const std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> key(hitset_key, hit_key);

  if (_do_cache)
  {
    const auto &  iter =
        _cache_get_energy_contribution_g4particle.find(make_pair(key, particle));
    if (iter != _cache_get_energy_contribution_g4particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hitset_key, hit_key);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
  {
    PHG4Hit* g4hit = *iter;
    if (get_truth_eval()->is_g4hit_from_particle(g4hit, particle))
    {
      energy += g4hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(
      make_pair(
          make_pair(
              key,
              particle),
          energy));

  return energy;
}

float SvtxHitEval::get_energy_contribution(const TrkrDefs::hitsetkey hitset_key, TrkrDefs::hitkey hit_key, PHG4Hit* g4hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }

  if (_strict)
  {
    assert(hit_key);
    assert(g4hit);
  }
  else if (!hit_key || !g4hit)
  {
    ++_errors;
    if (_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }

  const std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey> key(hitset_key, hit_key);
  if (_do_cache)
  {
    const auto &  iter =
        _cache_get_energy_contribution_g4hit.find(make_pair(key, g4hit));
    if (iter != _cache_get_energy_contribution_g4hit.end())
    {
      return iter->second;
    }
  }

  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.

  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hitset_key, hit_key);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
  {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(key, g4hit), energy));

  return energy;
}

void SvtxHitEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  for (const auto& [trkrID, trkrName] : TrkrDefs::TrkrNames)
  {
    TrkrHitSetContainer* _hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET_" + trkrName);
    if (_hitmap)
      _hitmaps[trkrID] = _hitmap;
  }

  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if (!_clustermap)
    _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");

  // need things off of the DST...
  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  _g4hits_mms = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return;
}

bool SvtxHitEval::has_node_pointers()
{
  if (_strict)
    assert(_hitmaps.size());
  else if (!_hitmaps.size())
  {
    return false;
  }

  if (_strict)
    assert(_g4hits_mms || _g4hits_tpc || _g4hits_intt || _g4hits_mvtx);
  else if (!_g4hits_mms && !_g4hits_tpc && !_g4hits_intt && !_g4hits_mvtx)
  {
    cout << "no hits" << endl;
    return false;
  }
  if (_strict)
    assert(_truthinfo);
  else if (!_truthinfo)
  {
    cout << " no truth" << endl;
    return false;
  }

  return true;
}
