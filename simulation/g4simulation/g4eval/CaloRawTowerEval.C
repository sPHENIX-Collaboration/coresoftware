
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
    _cache_all_towers_from_primary(),
    _cache_best_tower_from_primary(),
    _cache_get_energy_contribution_primary() {
}

void CaloRawTowerEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_primaries.clear();
  _cache_max_truth_primary_by_energy.clear();
  _cache_all_towers_from_primary.clear();
  _cache_best_tower_from_primary.clear();
  _cache_get_energy_contribution_primary.clear();

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
    cerr << PHWHERE << " ERROR: Can't find " << cellname << endl;
    exit(-1);
  }

  std::string hitname = "G4HIT_" + _caloname;  
  PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(_topNode,hitname.c_str());
  if (!g4hits){ 
    cerr << PHWHERE << " ERROR: Can't find " << hitname << endl;
    exit(-1);
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the towered cells
  std::pair< std::map<unsigned int,float>::const_iterator,
	     std::map<unsigned int,float>::const_iterator > cell_range = tower->get_g4cells();
  for (std::map<unsigned int, float>::const_iterator cell_iter = cell_range.first;
       cell_iter != cell_range.second; ++cell_iter) {
    unsigned int cell_id = cell_iter->first;
    PHG4CylinderCell *cell = g4cells->findCylinderCell(cell_id);

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
    PHG4Particle* primary = _trutheval.get_primary_particle( particle );
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
  std::set<PHG4Particle*> primaries = all_truth_primaries(tower);
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

std::set<RawTower*> CaloRawTowerEval::all_towers_from(PHG4Particle* primary) { 

  if (!_trutheval.is_primary(primary)) return std::set<RawTower*>();
  
  if ((_do_cache) &&
      (_cache_all_towers_from_primary.find(primary) != _cache_all_towers_from_primary.end())) {
    return _cache_all_towers_from_primary[primary];
  }
  
  // need things off of the DST...
  RawTowerContainer* towercontainer = findNode::getClass<RawTowerContainer>(_topNode,"RawTowerContainer");
  if (!towercontainer) {
    cerr << PHWHERE << " ERROR: Can't find RawTowerContainer" << endl;
    exit(-1);
  }

  std::set<RawTower*> towers;
  
  // loop over all the towers
  for (RawTowerContainer::Iterator iter = towercontainer->getTowers().first;
       iter != towercontainer->getTowers().second;
       ++iter) {

    RawTower* tower = iter->second;

    // loop over all truth particles connected to this tower
    std::set<PHG4Particle*> primaries = all_truth_primaries(tower);
    for (std::set<PHG4Particle*>::iterator jter = primaries.begin();
	 jter != primaries.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      if (candidate->get_track_id() == primary->get_track_id()) {
	towers.insert(tower);
      }    
    }
  }

  if (_do_cache) _cache_all_towers_from_primary.insert(make_pair(primary,towers));
  
  return towers;
}

RawTower* CaloRawTowerEval::best_tower_from(PHG4Particle* primary) {

  if (!_trutheval.is_primary(primary)) return NULL;
      
  if ((_do_cache) &&
      (_cache_best_tower_from_primary.find(primary) != _cache_best_tower_from_primary.end())) {
    return _cache_best_tower_from_primary[primary];
  }
  
  RawTower* best_tower = NULL;
  float best_energy = 0.0;  
  std::set<RawTower*> towers = all_towers_from(primary);
  for (std::set<RawTower*>::iterator iter = towers.begin();
       iter != towers.end();
       ++iter) {
    RawTower* tower = *iter;
    float energy = get_energy_contribution(tower,primary);
    if (energy > best_energy) {
      best_tower = tower;
      best_energy = energy;
    }
  }
 
  if (_do_cache) _cache_best_tower_from_primary.insert(make_pair(primary,best_tower));
  
  return best_tower;
}

// overlap calculations
float CaloRawTowerEval::get_energy_contribution(RawTower* tower, PHG4Particle* primary) {

  if (!_trutheval.is_primary(primary)) return NAN;
  
  if ((_do_cache) &&
      (_cache_get_energy_contribution_primary.find(make_pair(tower,primary)) !=
       _cache_get_energy_contribution_primary.end())) {
    return _cache_get_energy_contribution_primary[make_pair(tower,primary)];
  }
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(tower);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = *iter;
    if (g4hit->get_trkid() == primary->get_track_id()) {
      energy += g4hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_primary.insert(make_pair(make_pair(tower,primary),energy));
  
  return energy;
}
