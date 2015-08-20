
#include "CaloRawTowerEval.h"
#include "CaloTruthEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <string>
#include <cstdlib>
#include <set>
#include <map>
#include <float.h>

using namespace std;

CaloRawTowerEval::CaloRawTowerEval(PHCompositeNode* topNode, std::string caloname)
  : _topNode(topNode),
    _caloname(caloname),
    _trutheval(topNode,caloname),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_primaries(),
    _cache_max_truth_primary_by_energy(),
    _cache_all_hits_from_particle(),
    _cache_all_hits_from_g4hit(),
    _cache_best_hit_from_g4hit(),
    _cache_get_energy_contribution_g4particle(),
    _cache_get_energy_contribution_g4hit() {
}

void CaloRawTowerEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_max_truth_hit_by_energy.clear();
  _cache_all_truth_primaries.clear();
  _cache_max_truth_primary_by_energy.clear();
  _cache_all_hits_from_particle.clear();
  _cache_all_hits_from_g4hit.clear();
  _cache_best_hit_from_g4hit.clear();
  _cache_get_energy_contribution_g4particle.clear();
  _cache_get_energy_contribution_g4hit.clear();

  _trutheval.next_event(topNode);
  
  _topNode = topNode;  
}

std::set<PHG4Hit*> CaloRawTowerEval::all_truth_hits(RawTower* tower) {

  if ((_do_cache) &&
      (_cache_all_truth_hits.find(tower) != _cache_all_truth_hits.end())) {
    return _cache_all_truth_hits[tower];
  }
  
  // need things off of the DST...
  std::string cellname = "G4CELL_" + _caloname;
  PHG4CylinderCellContainer* g4cells = findNode::getClass<PHG4CylinderCellContainer>(_topNode,cellname.c_str());
  if (!g4cells) {
    cerr << PHWHERE << " ERROR: Can't find " << _cellname << endl;
    exit(-1);
  }

  std::string hitname = "G4HIT_" + _caloname;  
  PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(_topNode,hitname.c_str());
  if (!g4hits){ 
    cerr << PHWHERE << " ERROR: Can't find " << _hitname << endl;
    exit(-1);
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the towered cells
  std::pair< std::map<unsigned int,float>::const_iterator,
	     std::map<unsigned int,float>::const_iterator > cell_range = tower->get_g4cells();
  for (std::map<unsigned int, float>::const_iterator cell_iter = cell_range.first;
       cell_iter != cell_range.second; ++cell_iter) {
    unsigned int cell_id = cell_iter->first;
    PHG4CylinderCell *cell = g4cells->findCylinderCell(hit->get_cellid());

    // loop over all the g4hits in this cell
    for (PHG4CylinderCell::EdepConstIterator hit_iter = cell->get_g4hits().first;
	 hit_iter != cell->get_g4hits().second;
	 ++hit_iter) {      
      PHG4Hit* g4hit = g4hits->findHit(hit_iter->first);   
      // fill output set
      truth_hits.insert(g4hit);
    }
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(tower,truth_hits));
  
  return truth_hits;
}
  
std::set<PHG4Particle*> CaloRawTowerEval::all_truth_primaries(RawTower* tower) {

  if ((_do_cache) &&
      (_cache_all_truth_primaries.find(tower) != _cache_all_truth_primaries.end())) {
    return _cache_all_truth_primaries[tower];
  }
  
  std::set<PHG4Particle*> truth_primaries;
  
  std::set<PHG4Hit*> g4hits = all_truth_hits(tower);
  if (g4hits.empty()) return truth_primaries;
  
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(_topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = *iter;
    PHG4Particle* particle = truthinfo->GetHit( g4hit->get_trkid() );
    PHG4Particle* primary = _trutheval.get_primary(   
    if (!primary) continue;
    truth_primaries.insert(primary);
  }

  if (_do_cache) _cache_all_truth_primaries.insert(make_pair(tower,truth_primaries));
  
  return truth_primaries;
}

PHG4Particle* CaloRawTowerEval::max_truth_primary_by_energy(RawTower* tower) {

  if ((_do_cache) &&
      (_cache_max_truth_primary_by_energy.find(tower) != _cache_max_truth_primary_by_energy.end())) {
    return _cache_max_truth_primary_by_energy[tower];
  }
  
  // loop over all primaries associated with this tower and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_primary = NULL;
  float max_e = FLT_MIN;
  std::set<PHG4Particle*> primaries = all_truth_primaryies(tower);
  for (std::set<PHG4Particle*>::iterator iter = primaries.begin();
       iter != primaries.end();
       ++iter) {

    PHG4Particle* primary = *iter;
    float e = get_energy_contribution(tower,primary);
    if (e > max_e) {
      max_e = e;
      max_primary = primary;      
    }
  }

  if (_do_cache) _cache_max_truth_primary_by_energy.insert(make_pair(tower,max_primary));
  
  return max_primary;
}

std::set<RawTower*> CaloRawTowerEval::all_hits_from(PHG4Particle* g4particle) { 

  if ((_do_cache) &&
      (_cache_all_hits_from_particle.find(g4particle) != _cache_all_hits_from_particle.end())) {
    return _cache_all_hits_from_particle[g4particle];
  }
  
  // need things off of the DST...
  RawTowerMap* hitmap = findNode::getClass<RawTowerMap>(_topNode,"RawTowerMap");
  if (!hitmap) {
    cerr << PHWHERE << " ERROR: Can't find RawTowerMap" << endl;
    exit(-1);
  }

  std::set<RawTower*> hits;
  
  // loop over all the hits
  for (RawTowerMap::Iter iter = hitmap->begin();
       iter != hitmap->end();
       ++iter) {

    RawTower* hit = &iter->second;

    // loop over all truth particles connected to this hit
    std::set<PHG4Particle*> g4particles = all_truth_particles(hit);
    for (std::set<PHG4Particle*>::iterator jter = g4particles.begin();
	 jter != g4particles.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      if (candidate->get_track_id() == g4particle->get_track_id()) {
	hits.insert(hit);
      }    
    }
  }

  if (_do_cache) _cache_all_hits_from_particle.insert(make_pair(g4particle,hits));
  
  return hits;
}

std::set<RawTower*> CaloRawTowerEval::all_hits_from(PHG4Hit* g4hit) {

  if ((_do_cache) &&
      (_cache_all_hits_from_g4hit.find(g4hit) != _cache_all_hits_from_g4hit.end())) {
    return _cache_all_hits_from_g4hit[g4hit];
  }
  
  // need things off of the DST...
  RawTowerMap* hitmap = findNode::getClass<RawTowerMap>(_topNode,"RawTowerMap");
  if (!hitmap) {
    cerr << PHWHERE << " ERROR: Can't find RawTowerMap" << endl;
    exit(-1);
  }

  std::set<RawTower*> hits;
  
  // loop over all the hits
  for (RawTowerMap::Iter iter = hitmap->begin();
       iter != hitmap->end();
       ++iter) {

    RawTower* hit = &iter->second;

    // loop over all truth hits connected to this hit
    std::set<PHG4Hit*> g4hits = all_truth_hits(hit);
    for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
	 jter != g4hits.end();
	 ++jter) {
      PHG4Hit* candidate = *jter;
      if (candidate->get_hit_id() == g4hit->get_hit_id()) {
	hits.insert(hit);
      }    
    }
  }

  if (_do_cache) _cache_all_hits_from_g4hit.insert(make_pair(g4hit,hits));
  
  return hits;
}

RawTower* CaloRawTowerEval::best_hit_from(PHG4Hit* g4hit) {

  if ((_do_cache) &&
      (_cache_best_hit_from_g4hit.find(g4hit) != _cache_best_hit_from_g4hit.end())) {
    return _cache_best_hit_from_g4hit[g4hit];
  }
  
  // need things off of the DST...
  RawTowerMap* hitmap = findNode::getClass<RawTowerMap>(_topNode,"RawTowerMap");
  if (!hitmap) {
    cerr << PHWHERE << " ERROR: Can't find RawTowerMap" << endl;
    exit(-1);
  }

  RawTower* best_hit = NULL;
  float best_energy = 0.0;  
  std::set<RawTower*> hits = all_hits_from(g4hit);
  for (std::set<RawTower*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter) {
    RawTower* hit = *iter;
    float energy = get_energy_contribution(hit,g4hit);
    if (energy > best_energy) {
      best_hit = hit;
      best_energy = energy;
    }
  }
 
  if (_do_cache) _cache_best_hit_from_g4hit.insert(make_pair(g4hit,best_hit));
  
  return best_hit;
}

// overlap calculations
float CaloRawTowerEval::get_energy_contribution(RawTower* hit, PHG4Particle* particle) {

  if ((_do_cache) &&
      (_cache_get_energy_contribution_g4particle.find(make_pair(hit,particle)) !=
       _cache_get_energy_contribution_g4particle.end())) {
    return _cache_get_energy_contribution_g4particle[make_pair(hit,particle)];
  }
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = *iter;
    if (g4hit->get_trkid() == particle->get_track_id()) {
      energy += g4hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(make_pair(make_pair(hit,particle),energy));
  
  return energy;
}

float CaloRawTowerEval::get_energy_contribution(RawTower* hit, PHG4Hit* g4hit) {

  if ((_do_cache) &&
      (_cache_get_energy_contribution_g4hit.find(make_pair(hit,g4hit)) !=
       _cache_get_energy_contribution_g4hit.end())) {
    return _cache_get_energy_contribution_g4hit[make_pair(hit,g4hit)];
  }
  
  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;  
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(hit,g4hit),energy));
  
  return energy;
}
