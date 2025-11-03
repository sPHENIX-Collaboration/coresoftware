
#include "CaloRawTowerEval.h"
#include "CaloTruthEval.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/getClass.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>

class PHG4Hit;

CaloRawTowerEval::CaloRawTowerEval(PHCompositeNode* topNode, const std::string& caloname)
  : _caloname(caloname)
  , _trutheval(topNode, caloname)
{
  get_node_pointers(topNode);
}

CaloRawTowerEval::~CaloRawTowerEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      std::cout << "CaloRawTowerEval::~CaloRawTowerEval() - Error Count: " << _errors << std::endl;
    }
  }
}

void CaloRawTowerEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_primary_showers.clear();
  _cache_max_truth_primary_shower_by_energy.clear();
  _cache_all_towers_from_primary_shower.clear();
  _cache_best_tower_from_primary_shower.clear();
  _cache_get_energy_contribution_primary_shower.clear();

  _cache_all_truth_primary_particles.clear();
  _cache_max_truth_primary_particle_by_energy.clear();
  _cache_all_towers_from_primary_particle.clear();
  _cache_best_tower_from_primary_particle.clear();
  _cache_get_energy_contribution_primary_particle.clear();

  _cache_all_truth_hits.clear();

  _cache_towerinfo_all_truth_primary_showers.clear();
  _cache_towerinfo_max_truth_primary_shower_by_energy.clear();
  _cache_towerinfo_all_towers_from_primary_shower.clear();
  _cache_towerinfo_best_tower_from_primary_shower.clear();
  _cache_towerinfo_get_energy_contribution_primary_shower.clear();

  _cache_towerinfo_all_truth_primary_particles.clear();
  _cache_towerinfo_max_truth_primary_particle_by_energy.clear();
  _cache_towerinfo_all_towers_from_primary_particle.clear();
  _cache_towerinfo_best_tower_from_primary_particle.clear();
  _cache_towerinfo_get_energy_contribution_primary_particle.clear();

  _cache_towerinfo_all_truth_hits.clear();

  _trutheval.next_event(topNode);

  get_node_pointers(topNode);
}

bool CaloRawTowerEval::has_reduced_node_pointers()
{
  if (!get_truth_eval()->has_reduced_node_pointers())
  {
    return false;
  }

  if (_strict)
  {
    assert(_towers || _towerinfos);
  }
  else if (!_towers && !_towerinfos)
  {
    return false;
  }

  if (_strict)
  {
    assert(_truthinfo);
  }
  else if (!_truthinfo)
  {
    return false;
  }

  return true;
}

std::set<PHG4Shower*> CaloRawTowerEval::all_truth_primary_showers(RawTower* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_do_cache)
  {
    std::map<RawTower*, std::set<PHG4Shower*> >::iterator iter =
        _cache_all_truth_primary_showers.find(tower);
    if (iter != _cache_all_truth_primary_showers.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Shower*> showers;

  RawTower::ShowerConstRange shower_range = tower->get_g4showers();
  for (RawTower::ShowerConstIterator iter = shower_range.first;
       iter != shower_range.second;
       ++iter)
  {
    PHG4Shower* shower = _truthinfo->GetShower(iter->first);

    if (_strict)
    {
      assert(shower);
    }
    else if (!shower)
    {
      ++_errors;
      continue;
    }

    showers.insert(shower);
  }

  if (_do_cache)
  {
    _cache_all_truth_primary_showers.insert(std::make_pair(tower, showers));
  }

  return showers;
}

std::set<PHG4Shower*> CaloRawTowerEval::all_truth_primary_showers(TowerInfo* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_all_truth_primary_showers.find(tower); iter != _cache_towerinfo_all_truth_primary_showers.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Shower*> showers;
  // use const incase something bad happens
  const TowerInfo::ShowerEdepMap& showerEdepMap = tower->get_showerEdepMap();

  for (const auto& [showerID, edep] : showerEdepMap)
  {
    PHG4Shower* shower = _truthinfo->GetShower(showerID);

    // Check for shower validity
    if (_strict)
    {
      assert(shower);
    }
    else if (!shower)
    {
      ++_errors;
      continue;
    }

    // Process or store showers based on edep if needed
    showers.insert(shower);
  }

  if (_do_cache)
  {
    _cache_towerinfo_all_truth_primary_showers.insert(std::make_pair(tower, showers));
  }

  return showers;
}

PHG4Shower* CaloRawTowerEval::max_truth_primary_shower_by_energy(RawTower* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<RawTower*, PHG4Shower*>::iterator iter =
        _cache_max_truth_primary_shower_by_energy.find(tower);
    if (iter != _cache_max_truth_primary_shower_by_energy.end())
    {
      return iter->second;
    }
  }

  PHG4Shower* max_shower = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Shower*> showers = all_truth_primary_showers(tower);

  for (auto* shower : showers)
  {
    if (_strict)
    {
      assert(shower);
    }
    else if (!shower)
    {
      ++_errors;
      continue;
    }

    float e = get_energy_contribution(tower, shower);
    if (std::isnan(e))
    {
      continue;
    }
    if (e > max_e)
    {
      max_e = e;
      max_shower = shower;
    }
  }

  if (_do_cache)
  {
    _cache_max_truth_primary_shower_by_energy.insert(std::make_pair(tower, max_shower));
  }

  return max_shower;
}

PHG4Shower* CaloRawTowerEval::max_truth_primary_shower_by_energy(TowerInfo* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return nullptr;
  }
  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_max_truth_primary_shower_by_energy.find(tower); iter != _cache_towerinfo_max_truth_primary_shower_by_energy.end())
    {
      return iter->second;
    }
  }

  PHG4Shower* max_shower = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Shower*> showers = all_truth_primary_showers(tower);

  for (auto* shower : showers)
  {
    if (_strict)
    {
      assert(shower);
    }
    else if (!shower)
    {
      ++_errors;
      continue;
    }

    float e = get_energy_contribution(tower, shower);
    if (std::isnan(e))
    {
      continue;
    }
    if (e > max_e)
    {
      max_e = e;
      max_shower = shower;
    }
  }

  if (_do_cache)
  {
    _cache_towerinfo_max_truth_primary_shower_by_energy.insert(std::make_pair(tower, max_shower));
  }

  return max_shower;
}

RawTower* CaloRawTowerEval::best_tower_from(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++_errors;
    return nullptr;
  }

  if (!_trutheval.is_primary(shower))
  {
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Shower*, RawTower*>::iterator iter =
        _cache_best_tower_from_primary_shower.find(shower);
    if (iter != _cache_best_tower_from_primary_shower.end())
    {
      return iter->second;
    }
  }

  RawTower* best_tower = nullptr;
  float best_energy = FLT_MAX * -1.0;
  std::set<RawTower*> towers = all_towers_from(shower);
  for (auto* tower : towers)
  {
    if (_strict)
    {
      assert(tower);
    }
    else if (!tower)
    {
      ++_errors;
      continue;
    }

    float energy = get_energy_contribution(tower, shower);
    if (std::isnan(energy))
    {
      continue;
    }
    if (energy > best_energy)
    {
      best_tower = tower;
      best_energy = energy;
    }
  }

  if (_do_cache)
  {
    _cache_best_tower_from_primary_shower.insert(std::make_pair(shower, best_tower));
  }

  return best_tower;
}

TowerInfo* CaloRawTowerEval::best_towerinfo_from(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++_errors;
    return nullptr;
  }

  if (!_trutheval.is_primary(shower))
  {
    return nullptr;
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_best_tower_from_primary_shower.find(shower); iter != _cache_towerinfo_best_tower_from_primary_shower.end())
    {
      return iter->second;
    }
  }

  TowerInfo* best_tower = nullptr;
  float best_energy = FLT_MAX * -1.0;
  std::set<TowerInfo*> towers = all_towerinfos_from(shower);
  for (auto* tower : towers)
  {
    if (_strict)
    {
      assert(tower);
    }
    else if (!tower)
    {
      ++_errors;
      continue;
    }

    float energy = get_energy_contribution(tower, shower);
    if (std::isnan(energy))
    {
      continue;
    }
    if (energy > best_energy)
    {
      best_tower = tower;
      best_energy = energy;
    }
  }

  if (_do_cache)
  {
    _cache_towerinfo_best_tower_from_primary_shower.insert(std::make_pair(shower, best_tower));
  }

  return best_tower;
}

std::set<RawTower*> CaloRawTowerEval::all_towers_from(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<RawTower*>();
  }

  if (_strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++_errors;
    return std::set<RawTower*>();
  }

  if (!_trutheval.is_primary(shower))
  {
    return std::set<RawTower*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Shower*, std::set<RawTower*> >::iterator iter =
        _cache_all_towers_from_primary_shower.find(shower);
    if (iter != _cache_all_towers_from_primary_shower.end())
    {
      return iter->second;
    }
  }

  std::set<RawTower*> towers;

  // loop over all the towers
  for (RawTowerContainer::Iterator iter = _towers->getTowers().first;
       iter != _towers->getTowers().second;
       ++iter)
  {
    RawTower* tower = iter->second;

    std::set<PHG4Shower*> showers = all_truth_primary_showers(tower);
    for (auto* candidate : showers)
    {
      if (_strict)
      {
        assert(candidate);
      }
      else if (!candidate)
      {
        ++_errors;
        continue;
      }

      if (candidate->get_id() == shower->get_id())
      {
        towers.insert(tower);
      }
    }
  }

  if (_do_cache)
  {
    _cache_all_towers_from_primary_shower.insert(std::make_pair(shower, towers));
  }

  return towers;
}

std::set<TowerInfo*> CaloRawTowerEval::all_towerinfos_from(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<TowerInfo*>();
  }

  if (_strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++_errors;
    return std::set<TowerInfo*>();
  }

  if (!_trutheval.is_primary(shower))
  {
    return std::set<TowerInfo*>();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_all_towers_from_primary_shower.find(shower); iter != _cache_towerinfo_all_towers_from_primary_shower.end())
    {
      return iter->second;
    }
  }

  std::set<TowerInfo*> towers;

  // loop over all the towers
  unsigned int ntowers = _towerinfos->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo* tower = _towerinfos->get_tower_at_channel(channel);

    std::set<PHG4Shower*> showers = all_truth_primary_showers(tower);
    for (auto* candidate : showers)
    {
      if (_strict)
      {
        assert(candidate);
      }
      else if (!candidate)
      {
        ++_errors;
        continue;
      }

      if (candidate->get_id() == shower->get_id())
      {
        towers.insert(tower);
      }
    }
  }

  if (_do_cache)
  {
    _cache_towerinfo_all_towers_from_primary_shower.insert(std::make_pair(shower, towers));
  }

  return towers;
}

float CaloRawTowerEval::get_energy_contribution(RawTower* tower, PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_strict)
  {
    assert(tower);
    assert(shower);
  }
  else if (!tower || !shower)
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (!_trutheval.is_primary(shower))
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_do_cache)
  {
    std::map<std::pair<RawTower*, PHG4Shower*>, float>::iterator iter =
        _cache_get_energy_contribution_primary_shower.find(std::make_pair(tower, shower));
    if (iter != _cache_get_energy_contribution_primary_shower.end())
    {
      return iter->second;
    }
  }

  // loop over the tower shower entries
  float energy = 0.0;
  RawTower::ShowerConstRange range = tower->get_g4showers();
  RawTower::ShowerConstIterator iter = tower->find_g4shower(shower->get_id());
  if (iter != range.second)
  {
    energy = iter->second;
  }

  if (_do_cache)
  {
    _cache_get_energy_contribution_primary_shower.insert(std::make_pair(std::make_pair(tower, shower), energy));
  }

  return energy;
}

float CaloRawTowerEval::get_energy_contribution(TowerInfo* tower, PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_strict)
  {
    assert(tower);
    assert(shower);
  }
  else if (!tower || !shower)
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (!_trutheval.is_primary(shower))
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_get_energy_contribution_primary_shower.find(std::make_pair(tower, shower)); iter != _cache_towerinfo_get_energy_contribution_primary_shower.end())
    {
      return iter->second;
    }
  }

  // loop over the tower shower entries
  float energy = 0.0;

  const TowerInfo::ShowerEdepMap& showerEdepMap = tower->get_showerEdepMap();
  auto showerIter = showerEdepMap.find(shower->get_id());
  if (showerIter != showerEdepMap.end())
  {
    energy = showerIter->second;
  }

  if (_do_cache)
  {
    _cache_towerinfo_get_energy_contribution_primary_shower.insert(std::make_pair(std::make_pair(tower, shower), energy));
  }

  return energy;
}

std::set<PHG4Particle*> CaloRawTowerEval::all_truth_primary_particles(RawTower* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<RawTower*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_primary_particles.find(tower);
    if (iter != _cache_all_truth_primary_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_primaries;

  std::set<PHG4Shower*> showers = all_truth_primary_showers(tower);

  for (auto* shower : showers)
  {
    PHG4Particle* primary = get_truth_eval()->get_primary_particle(shower);

    if (_strict)
    {
      assert(primary);
    }
    else if (!primary)
    {
      ++_errors;
      continue;
    }

    truth_primaries.insert(primary);
  }

  if (_do_cache)
  {
    _cache_all_truth_primary_particles.insert(std::make_pair(tower, truth_primaries));
  }

  return truth_primaries;
}

std::set<PHG4Particle*> CaloRawTowerEval::all_truth_primary_particles(TowerInfo* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_all_truth_primary_particles.find(tower); iter != _cache_towerinfo_all_truth_primary_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_primaries;

  std::set<PHG4Shower*> showers = all_truth_primary_showers(tower);

  for (auto* shower : showers)
  {
    PHG4Particle* primary = get_truth_eval()->get_primary_particle(shower);

    if (_strict)
    {
      assert(primary);
    }
    else if (!primary)
    {
      ++_errors;
      continue;
    }

    truth_primaries.insert(primary);
  }

  if (_do_cache)
  {
    _cache_towerinfo_all_truth_primary_particles.insert(std::make_pair(tower, truth_primaries));
  }

  return truth_primaries;
}

PHG4Particle* CaloRawTowerEval::max_truth_primary_particle_by_energy(RawTower* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<RawTower*, PHG4Particle*>::iterator iter =
        _cache_max_truth_primary_particle_by_energy.find(tower);
    if (iter != _cache_max_truth_primary_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  PHG4Particle* max_primary = nullptr;
  PHG4Shower* max_shower = max_truth_primary_shower_by_energy(tower);

  if (max_shower)
  {
    max_primary = get_truth_eval()->get_primary_particle(max_shower);
  }

  if (_do_cache)
  {
    _cache_max_truth_primary_particle_by_energy.insert(std::make_pair(tower, max_primary));
  }

  return max_primary;
}

PHG4Particle* CaloRawTowerEval::max_truth_primary_particle_by_energy(TowerInfo* tower)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_max_truth_primary_particle_by_energy.find(tower); iter != _cache_towerinfo_max_truth_primary_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  PHG4Particle* max_primary = nullptr;
  PHG4Shower* max_shower = max_truth_primary_shower_by_energy(tower);

  if (max_shower)
  {
    max_primary = get_truth_eval()->get_primary_particle(max_shower);
  }

  if (_do_cache)
  {
    _cache_towerinfo_max_truth_primary_particle_by_energy.insert(std::make_pair(tower, max_primary));
  }

  return max_primary;
}

std::set<RawTower*> CaloRawTowerEval::all_towers_from(PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<RawTower*>();
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<RawTower*>();
  }

  if (!_trutheval.is_primary(primary))
  {
    return std::set<RawTower*>();
  }

  // use primary map pointer
  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<RawTower*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<RawTower*> >::iterator iter =
        _cache_all_towers_from_primary_particle.find(primary);
    if (iter != _cache_all_towers_from_primary_particle.end())
    {
      return iter->second;
    }
  }

  std::set<RawTower*> towers;

  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);

  if (shower)
  {
    towers = all_towers_from(shower);
  }

  if (_do_cache)
  {
    _cache_all_towers_from_primary_particle.insert(std::make_pair(primary, towers));
  }

  return towers;
}

std::set<TowerInfo*> CaloRawTowerEval::all_towerinfos_from(PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<TowerInfo*>();
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<TowerInfo*>();
  }

  if (!_trutheval.is_primary(primary))
  {
    return std::set<TowerInfo*>();
  }

  // use primary map pointer
  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<TowerInfo*>();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_all_towers_from_primary_particle.find(primary); iter != _cache_towerinfo_all_towers_from_primary_particle.end())
    {
      return iter->second;
    }
  }

  std::set<TowerInfo*> towers;

  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);

  if (shower)
  {
    towers = all_towerinfos_from(shower);
  }

  if (_do_cache)
  {
    _cache_towerinfo_all_towers_from_primary_particle.insert(std::make_pair(primary, towers));
  }

  return towers;
}

RawTower* CaloRawTowerEval::best_tower_from(PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return nullptr;
  }

  if (!_trutheval.is_primary(primary))
  {
    return nullptr;
  }

  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, RawTower*>::iterator iter =
        _cache_best_tower_from_primary_particle.find(primary);
    if (iter != _cache_best_tower_from_primary_particle.end())
    {
      return iter->second;
    }
  }

  RawTower* best_tower = nullptr;
  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);
  if (shower)
  {
    best_tower = best_tower_from(shower);
  }

  if (_do_cache)
  {
    _cache_best_tower_from_primary_particle.insert(std::make_pair(primary, best_tower));
  }

  return best_tower;
}

TowerInfo* CaloRawTowerEval::best_towerinfo_from(PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return nullptr;
  }

  if (!_trutheval.is_primary(primary))
  {
    return nullptr;
  }

  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_best_tower_from_primary_particle.find(primary); iter != _cache_towerinfo_best_tower_from_primary_particle.end())
    {
      return iter->second;
    }
  }

  TowerInfo* best_tower = nullptr;
  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);
  if (shower)
  {
    best_tower = best_towerinfo_from(shower);
  }

  if (_do_cache)
  {
    _cache_towerinfo_best_tower_from_primary_particle.insert(std::make_pair(primary, best_tower));
  }

  return best_tower;
}

// overlap calculations
float CaloRawTowerEval::get_energy_contribution(RawTower* tower, PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_strict)
  {
    assert(tower);
    assert(primary);
  }
  else if (!tower || !primary)
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (!_trutheval.is_primary(primary))
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  // reduce cache misses by using only pointer from PrimaryMap
  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_do_cache)
  {
    std::map<std::pair<RawTower*, PHG4Particle*>, float>::iterator iter =
        _cache_get_energy_contribution_primary_particle.find(std::make_pair(tower, primary));
    if (iter != _cache_get_energy_contribution_primary_particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;

  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);

  if (shower)
  {
    energy = get_energy_contribution(tower, shower);
  }

  if (_do_cache)
  {
    _cache_get_energy_contribution_primary_particle.insert(std::make_pair(std::make_pair(tower, primary), energy));
  }

  return energy;
}

float CaloRawTowerEval::get_energy_contribution(TowerInfo* tower, PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_strict)
  {
    assert(tower);
    assert(primary);
  }
  else if (!tower || !primary)
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (!_trutheval.is_primary(primary))
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  // reduce cache misses by using only pointer from PrimaryMap
  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::numeric_limits<float>::quiet_NaN();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_get_energy_contribution_primary_particle.find(std::make_pair(tower, primary)); iter != _cache_towerinfo_get_energy_contribution_primary_particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;

  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);

  if (shower)
  {
    energy = get_energy_contribution(tower, shower);
  }

  if (_do_cache)
  {
    _cache_towerinfo_get_energy_contribution_primary_particle.insert(std::make_pair(std::make_pair(tower, primary), energy));
  }

  return energy;
}

bool CaloRawTowerEval::has_full_node_pointers()
{
  if (!get_truth_eval()->has_full_node_pointers())
  {
    return false;
  }

  if (_strict)
  {
    assert(_towers || _towerinfos);
  }
  else if (!_towers && !_towerinfos)
  {
    return false;
  }

  if (_strict)
  {
    assert(_g4cells);
  }
  else if (!_g4cells)
  {
    return false;
  }

  if (_strict)
  {
    assert(_g4hits);
  }
  else if (!_g4hits)
  {
    return false;
  }

  if (_strict)
  {
    assert(_truthinfo);
  }
  else if (!_truthinfo)
  {
    return false;
  }

  return true;
}

std::set<PHG4Hit*> CaloRawTowerEval::all_truth_hits(RawTower* tower)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<RawTower*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(tower);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;
  // loop over all the towered cells
  RawTower::CellConstRange cell_range = tower->get_g4cells();
  for (RawTower::CellConstIterator cell_iter = cell_range.first;
       cell_iter != cell_range.second; ++cell_iter)
  {
    unsigned int cell_id = cell_iter->first;
    PHG4Cell* cell = _g4cells->findCell(cell_id);

    if (_strict)
    {
      assert(cell);
    }
    else if (!cell)
    {
      ++_errors;
      continue;
    }

    // loop over all the g4hits in this cell
    for (PHG4Cell::EdepConstIterator hit_iter = cell->get_g4hits().first;
         hit_iter != cell->get_g4hits().second;
         ++hit_iter)
    {
      PHG4Hit* g4hit = _g4hits->findHit(hit_iter->first);

      if (_strict)
      {
        assert(g4hit);
      }
      else if (!g4hit)
      {
        ++_errors;
        continue;
      }

      // fill output set
      truth_hits.insert(g4hit);
    }
  }

  if (_do_cache)
  {
    _cache_all_truth_hits.insert(std::make_pair(tower, truth_hits));
  }

  return truth_hits;
}

std::set<PHG4Hit*> CaloRawTowerEval::all_truth_hits(TowerInfo* tower)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(tower);
  }
  else if (!tower)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    if (auto iter = _cache_towerinfo_all_truth_hits.find(tower); iter != _cache_towerinfo_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;
  const TowerInfo::EdepMap& edepMap = tower->get_hitEdepMap();
  // loop over all g4hits
  for (const auto& [hitID, edep] : edepMap)
  {
    PHG4Hit* g4hit = _g4hits->findHit(hitID);

    if (_strict)
    {
      assert(g4hit);
    }
    else if (!g4hit)
    {
      ++_errors;
      continue;
    }
    // fill output set
    truth_hits.insert(g4hit);
  }

  if (_do_cache)
  {
    _cache_towerinfo_all_truth_hits.insert(std::make_pair(tower, truth_hits));
  }

  return truth_hits;
}

void CaloRawTowerEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  std::string towername = "TOWER_CALIB_" + _caloname;
  _towers = findNode::getClass<RawTowerContainer>(topNode, towername);

  std::string towerinfoname = "TOWERINFO_CALIB_" + _caloname;
  _towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerinfoname);

  std::string cellname = "G4CELL_" + _caloname;
  _g4cells = findNode::getClass<PHG4CellContainer>(topNode, cellname);

  std::string hitname = "G4HIT_" + _caloname;
  _g4hits = findNode::getClass<PHG4HitContainer>(topNode, hitname);

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return;
}
