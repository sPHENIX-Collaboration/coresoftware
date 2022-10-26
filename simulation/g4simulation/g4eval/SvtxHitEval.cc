#include "SvtxHitEval.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>
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
#include <iostream>                          // for operator<<, endl, basic_...
#include <map>
#include <set>

class TrkrHit;

using namespace std;

SvtxHitEval::SvtxHitEval(PHCompositeNode* topNode)
  : _trutheval(topNode)
  , _hitmap(nullptr)
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

/*
PHG4Cell* SvtxHitEval::get_cell(SvtxHit* hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit);
  }
  else if (!hit)
  {
    ++_errors;
    cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  // hop from reco hit to g4cell
  PHG4Cell* cell = nullptr;
  if (_g4cells_svtx) cell = _g4cells_svtx->findCell(hit->get_cellid());
  if (!cell && _g4cells_tracker) cell = _g4cells_tracker->findCell(hit->get_cellid());
  if (!cell && _g4cells_maps) cell = _g4cells_maps->findCell(hit->get_cellid());

  // only noise hits (cellid left at default value) should not trace
  if ((_strict) && (hit->get_cellid() != 0xFFFFFFFF))
  {
    assert(cell);
  }
  else if (!cell)
  {
    ++_errors;
    cout << PHWHERE << " nerr: " << _errors << endl;
  }

  return cell;
}
*/

std::set<PHG4Hit*> SvtxHitEval::all_truth_hits(TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(hit_key);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  /*
  // hop from reco hit to g4cell
  PHG4Cell* cell = nullptr;
  if (_g4cells_svtx) cell = _g4cells_svtx->findCell(hit->get_cellid());
  if (!cell && _g4cells_tracker) cell = _g4cells_tracker->findCell(hit->get_cellid());
  if (!cell && _g4cells_maps) cell = _g4cells_maps->findCell(hit->get_cellid());

  // only noise hits (cellid left at default value) should not trace
  if ((_strict) && (hit->get_cellid() != 0xFFFFFFFF))
    assert(cell);
  else if (!cell)
  {
    ++_errors;
    cout << PHWHERE << " nerr: " << _errors << endl;
    return truth_hits;
  }

  //cout << "Eval: hitid " << hit->get_id() << " cellid " << cell->get_cellid() << endl;
  // loop over all the g4hits in this cell
  for (PHG4Cell::EdepConstIterator g4iter = cell->get_g4hits().first;
       g4iter != cell->get_g4hits().second;
       ++g4iter)
  {
    //cout << "    Looking for hit " << g4iter->first << " in layer " << cell->get_layer() << " with edep " << g4iter->second << endl;
    PHG4Hit* g4hit = nullptr;
    if (_g4hits_svtx) g4hit = _g4hits_svtx->findHit(g4iter->first);
    if (!g4hit && _g4hits_tracker) g4hit = _g4hits_tracker->findHit(g4iter->first);
    if (!g4hit && _g4hits_maps) g4hit = _g4hits_maps->findHit(g4iter->first);
    if (!g4hit) cout << "    Failed to find  g4hit " << g4iter->first << " with edep " << g4iter->second << endl;
    if (_strict)
      assert(g4hit);
    else if (!g4hit)
    {
      ++_errors;
      cout << PHWHERE << " nerr: " << _errors << endl;
      continue;
    }
  */

  // get all of the g4hits for this hit_key
  // have to start with all hitsets, unfortunately
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets(); 
  for(TrkrHitSetContainer::ConstIterator iter = all_hitsets.first; iter != all_hitsets.second; ++iter)
    {
      TrkrDefs::hitsetkey hitset_key = iter->first;
      unsigned int trkrid = TrkrDefs::getTrkrId(hitset_key);
      TrkrHitSet *hitset = iter->second;

      // does this hitset contain our hitkey?
      TrkrHit *hit = nullptr;
      hit = hitset->getHit(hit_key);
      if(hit)
	{
	  // get g4hits for this hit

	  std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
	  _hit_truth_map->getG4Hits(hitset_key, hit_key, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
	  for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); 
	      htiter != temp_map.end(); ++htiter) 
	    {
	      // extract the g4 hit key here and add the g4hit to the set
	      PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	      //cout << "           hitkey " << hitkey <<  " g4hitkey " << g4hitkey << endl;	  
	      PHG4Hit * g4hit = nullptr;
       switch( trkrid )
       {
        case TrkrDefs::tpcId: g4hit = _g4hits_tpc->findHit(g4hitkey); break;
        case TrkrDefs::inttId: g4hit = _g4hits_intt->findHit(g4hitkey); break;
        case TrkrDefs::mvtxId: g4hit = _g4hits_mvtx->findHit(g4hitkey); break;
       case TrkrDefs::micromegasId: g4hit = _g4hits_mms->findHit(g4hitkey); break;
        default: break;
       }
	      // fill output set
	      if( g4hit ) truth_hits.insert(g4hit);
	    }
	}
    }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(hit_key, truth_hits));

  return truth_hits;
}

std::set<PHG4Hit*> SvtxHitEval::all_truth_hits(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(hit_key);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets(trkrid); 

  for(TrkrHitSetContainer::ConstIterator iter = all_hitsets.first; iter != all_hitsets.second; ++iter)
  {
    TrkrDefs::hitsetkey hitset_key = iter->first;
    TrkrHitSet *hitset = iter->second;

    // does this hitset contain our hitkey?
    TrkrHit *hit = nullptr;
    hit = hitset->getHit(hit_key);
    if(hit)
      {
        // get g4hits for this hit

        std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
        _hit_truth_map->getG4Hits(hitset_key, hit_key, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
        for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); 
            htiter != temp_map.end(); ++htiter) 
          {
            // extract the g4 hit key here and add the g4hit to the set
            PHG4HitDefs::keytype g4hitkey = htiter->second.second;
            //cout << "           hitkey " << hitkey <<  " g4hitkey " << g4hitkey << endl;	  
            PHG4Hit * g4hit = nullptr;
     switch( trkrid )
     {
      case TrkrDefs::tpcId: g4hit = _g4hits_tpc->findHit(g4hitkey); break;
      case TrkrDefs::inttId: g4hit = _g4hits_intt->findHit(g4hitkey); break;
      case TrkrDefs::mvtxId: g4hit = _g4hits_mvtx->findHit(g4hitkey); break;
     case TrkrDefs::micromegasId: g4hit = _g4hits_mms->findHit(g4hitkey); break;
      default: break;
     }
            // fill output set
            if( g4hit ) truth_hits.insert(g4hit);
          }
      }
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(hit_key, truth_hits));

  return truth_hits;
}

PHG4Hit* SvtxHitEval::max_truth_hit_by_energy(TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, PHG4Hit*>::iterator iter =
        _cache_max_truth_hit_by_energy.find(hit_key);
    if (iter != _cache_max_truth_hit_by_energy.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> hits = all_truth_hits(hit_key);
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

  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(hit_key, max_hit));

  return max_hit;
}

PHG4Hit* SvtxHitEval::max_truth_hit_by_energy(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, PHG4Hit*>::iterator iter =
        _cache_max_truth_hit_by_energy.find(hit_key);
    if (iter != _cache_max_truth_hit_by_energy.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> hits = all_truth_hits(hit_key, trkrid);
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

  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(hit_key, max_hit));

  return max_hit;
}

std::set<PHG4Particle*> SvtxHitEval::all_truth_particles(TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(hit_key);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key);

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
      if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
      continue;
    }

    truth_particles.insert(particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(hit_key, truth_particles));

  return truth_particles;
}

std::set<PHG4Particle*> SvtxHitEval::all_truth_particles(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0) cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(hit_key);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key, trkrid);

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
      if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
      continue;
    }

    truth_particles.insert(particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(hit_key, truth_particles));

  return truth_particles;
}

PHG4Particle* SvtxHitEval::max_truth_particle_by_energy(TrkrDefs::hitkey hit_key)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, PHG4Particle*>::iterator iter =
        _cache_max_truth_particle_by_energy.find(hit_key);
    if (iter != _cache_max_truth_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  // loop over all particles associated with this hit and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(hit_key);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(hit_key, particle);
    if (e > max_e)
    {
      max_e = e;
      max_particle = particle;
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(hit_key, max_particle));

  return max_particle;
}

PHG4Particle* SvtxHitEval::max_truth_particle_by_energy(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_strict)
  {
    assert(hit_key);
  }
  else if (!hit_key)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<TrkrDefs::hitkey, PHG4Particle*>::iterator iter =
        _cache_max_truth_particle_by_energy.find(hit_key);
    if (iter != _cache_max_truth_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  // loop over all particles associated with this hit and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(hit_key, trkrid);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(hit_key, particle);
    if (e > max_e)
    {
      max_e = e;
      max_particle = particle;
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(hit_key, max_particle));

  return max_particle;
}

std::set<TrkrDefs::hitkey> SvtxHitEval::all_hits_from(PHG4Particle* g4particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
     if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<TrkrDefs::hitkey>();
  }

  if (_strict)
  {
    assert(g4particle);
  }
  else if (!g4particle)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<TrkrDefs::hitkey>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<TrkrDefs::hitkey> >::iterator iter =
        _cache_all_hits_from_particle.find(g4particle);
    if (iter != _cache_all_hits_from_particle.end())
    {
      return iter->second;
    }
  }

  std::set<TrkrDefs::hitkey> hits;

  // loop over all the hits
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
       iter != all_hitsets.second;
       ++iter)
  {    
    TrkrHitSet::ConstRange range = iter->second->getHits();
    for(TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
      {
	TrkrDefs::hitkey hit_key = hitr->first;
	
	// loop over all truth hits connected to this hit
	std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key);
	for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
	     jter != g4hits.end();
	     ++jter)
	  {
	    PHG4Hit* candidate = *jter;
	    PHG4Particle *particle = _truthinfo->GetParticle(candidate->get_trkid());
	    if (g4particle->get_track_id() == particle->get_track_id()) 
	      {
		hits.insert(hit_key);
	      }
	  }
      }	
  }
  
  if (_do_cache) _cache_all_hits_from_particle.insert(make_pair(g4particle, hits));
  
  return hits;
}

std::set<TrkrDefs::hitkey> SvtxHitEval::all_hits_from(PHG4Hit* g4hit)
{
  if (!has_node_pointers())
    {
    ++_errors;
     if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<TrkrDefs::hitkey>();
  }

  if (_strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return std::set<TrkrDefs::hitkey>();
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, std::set<TrkrDefs::hitkey> >::iterator iter =
        _cache_all_hits_from_g4hit.find(g4hit);
    if (iter != _cache_all_hits_from_g4hit.end())
    {
      return iter->second;
    }
  }

  std::set<TrkrDefs::hitkey> hits;

  unsigned int hit_layer = g4hit->get_layer();

  // loop over all the hits
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter = all_hitsets.first;
       iter != all_hitsets.second;
       ++iter)
  {
    TrkrHitSet::ConstRange range = iter->second->getHits();
    for(TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
      {
	TrkrDefs::hitkey hit_key = hitr->first;

	if (TrkrDefs::getLayer(hit_key) != hit_layer) continue;

	// loop over all truth hits connected to this hit
	std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key);
	for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
	     jter != g4hits.end();
	     ++jter)
	  {
	    PHG4Hit* candidate = *jter;
	    if (candidate->get_hit_id() == g4hit->get_hit_id())
	      {
		hits.insert(hit_key);
	      }
	  }
      }
  }
  
  if (_do_cache) _cache_all_hits_from_g4hit.insert(make_pair(g4hit, hits));
  
  return hits;
}

 TrkrDefs::hitkey SvtxHitEval::best_hit_from(PHG4Hit* g4hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)   cout << PHWHERE << " nerr: " << _errors << endl;
    return 0;
  }

  if (_strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return 0;
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, TrkrDefs::hitkey>::iterator iter =
        _cache_best_hit_from_g4hit.find(g4hit);
    if (iter != _cache_best_hit_from_g4hit.end())
    {
      return iter->second;
    }
  }

  TrkrDefs::hitkey best_hit = 0;
  float best_energy = 0.0;
  std::set<TrkrDefs::hitkey> hits = all_hits_from(g4hit);
  for (std::set<TrkrDefs::hitkey>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter)
  {
    TrkrDefs::hitkey hit_key = *iter;
    float energy = get_energy_contribution(hit_key, g4hit);
    if (energy > best_energy)
    {
      best_hit = hit_key;
      best_energy = energy;
    }
  }

  if (_do_cache) _cache_best_hit_from_g4hit.insert(make_pair(g4hit, best_hit));

  return best_hit;
}

// overlap calculations
 float SvtxHitEval::get_energy_contribution(TrkrDefs::hitkey hit_key, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
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
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<std::pair<TrkrDefs::hitkey, PHG4Particle*>, float>::iterator iter =
        _cache_get_energy_contribution_g4particle.find(make_pair(hit_key, particle));
    if (iter != _cache_get_energy_contribution_g4particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key);
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

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(make_pair(make_pair(hit_key, particle), energy));

  return energy;
}

 float SvtxHitEval::get_energy_contribution(TrkrDefs::hitkey hit_key, PHG4Hit* g4hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
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
    if(_verbosity > 0)  cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<std::pair<TrkrDefs::hitkey, PHG4Hit*>, float>::iterator iter =
        _cache_get_energy_contribution_g4hit.find(make_pair(hit_key, g4hit));
    if (iter != _cache_get_energy_contribution_g4hit.end())
    {
      return iter->second;
    }
  }

  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.

  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit_key);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
  {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(hit_key, g4hit), energy));

  return energy;
}

void SvtxHitEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  _hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
 
  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if(!_clustermap)
    _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  
 _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

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
    assert(_hitmap);
  else if (!_hitmap)
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
