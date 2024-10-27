
#include "CaloRawClusterEval.h"
#include "CaloTruthEval.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <phool/getClass.h>

#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>

class RawTower;

CaloRawClusterEval::CaloRawClusterEval(PHCompositeNode* topNode, const std::string& caloname)
  : _caloname(caloname)
  , _towereval(topNode, caloname)
{
  get_node_pointers(topNode);
}

CaloRawClusterEval::~CaloRawClusterEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      std::cout << "CaloRawClusterEval::~CaloRawClusterEval() - Error Count: " << _errors << std::endl;
    }
  }
}

void CaloRawClusterEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_primary_showers.clear();
  _cache_max_truth_primary_shower_by_energy.clear();
  _cache_all_clusters_from_primary_shower.clear();
  _cache_best_cluster_from_primary_shower.clear();
  _cache_get_energy_contribution_primary_shower.clear();

  _cache_all_truth_primary_particles.clear();
  _cache_max_truth_primary_particle_by_energy.clear();
  _cache_all_clusters_from_primary_particle.clear();
  _cache_best_cluster_from_primary_particle.clear();
  _cache_get_energy_contribution_primary_particle.clear();

  _cache_all_truth_hits.clear();

  _towereval.next_event(topNode);

  get_node_pointers(topNode);
}

bool CaloRawClusterEval::has_reduced_node_pointers()
{
  if (!get_rawtower_eval()->has_reduced_node_pointers())
  {
    return false;
  }

  if (_strict)
  {
    assert(_clusters);
  }
  else if (!_clusters)
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

  return true;
}

std::set<PHG4Shower*> CaloRawClusterEval::all_truth_primary_showers(RawCluster* cluster)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_do_cache)
  {
    std::map<RawCluster*, std::set<PHG4Shower*> >::iterator iter =
        _cache_all_truth_primary_showers.find(cluster);
    if (iter != _cache_all_truth_primary_showers.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Shower*> truth_primary_showers;

  // loop over all the clustered towers
  if (_usetowerinfo)
  {
    const RawCluster::TowerMap& tower_map = cluster->get_towermap();
    for (auto tower_iter : tower_map)
    {
      RawTowerDefs::keytype tower_key = tower_iter.first;
      unsigned int towerinfo_key = get_towerinfo_key(tower_key);
      TowerInfo* tower = _towerinfos->get_tower_at_key(towerinfo_key);
      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      std::set<PHG4Shower*> new_primary_showers = _towereval.all_truth_primary_showers(tower);

      for (auto shower : new_primary_showers)
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

        truth_primary_showers.insert(shower);
      }
    }
  }
  else
  {
    RawCluster::TowerConstRange begin_end = cluster->get_towers();
    for (RawCluster::TowerConstIterator iter = begin_end.first;
         iter != begin_end.second;
         ++iter)
    {
      RawTower* tower = _towers->getTower(iter->first);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      std::set<PHG4Shower*> new_primary_showers = _towereval.all_truth_primary_showers(tower);

      for (auto shower : new_primary_showers)
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

        truth_primary_showers.insert(shower);
      }
    }
  }
  if (_do_cache)
  {
    _cache_all_truth_primary_showers.insert(std::make_pair(cluster, truth_primary_showers));
  }

  return truth_primary_showers;
}

PHG4Shower* CaloRawClusterEval::max_truth_primary_shower_by_energy(RawCluster* cluster)
{
  if (!has_reduced_node_pointers())
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
    std::map<RawCluster*, PHG4Shower*>::iterator iter =
        _cache_max_truth_primary_shower_by_energy.find(cluster);
    if (iter != _cache_max_truth_primary_shower_by_energy.end())
    {
      return iter->second;
    }
  }

  // loop over all primaries associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Shower* max_primary = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Shower*> primary_showers = all_truth_primary_showers(cluster);
  for (auto primary : primary_showers)
  {
    if (_strict)
    {
      assert(primary);
    }
    else if (!primary)
    {
      ++_errors;
      continue;
    }

    float e = get_energy_contribution(cluster, primary);
    if (isnan(e))
    {
      continue;
    }
    if (e > max_e)
    {
      max_e = e;
      max_primary = primary;
    }
  }

  if (_do_cache)
  {
    _cache_max_truth_primary_shower_by_energy.insert(std::make_pair(cluster, max_primary));
  }

  return max_primary;
}

std::set<RawCluster*> CaloRawClusterEval::all_clusters_from(PHG4Shower* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<RawCluster*>();
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<RawCluster*>();
  }

  if (!get_truth_eval()->is_primary(primary))
  {
    return std::set<RawCluster*>();
  }

  primary = get_truth_eval()->get_primary_shower(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<RawCluster*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Shower*, std::set<RawCluster*> >::iterator iter =
        _cache_all_clusters_from_primary_shower.find(primary);
    if (iter != _cache_all_clusters_from_primary_shower.end())
    {
      return iter->second;
    }
  }

  std::set<RawCluster*> clusters;

  // loop over all the clusters
  for (RawClusterContainer::Iterator iter = _clusters->getClusters().first;
       iter != _clusters->getClusters().second;
       ++iter)
  {
    RawCluster* cluster = iter->second;

    std::set<PHG4Shower*> primary_showers = all_truth_primary_showers(cluster);
    for (auto candidate : primary_showers)
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

      if (get_truth_eval()->are_same_shower(candidate, primary))
      {
        clusters.insert(cluster);
      }
    }
  }

  if (_do_cache)
  {
    _cache_all_clusters_from_primary_shower.insert(std::make_pair(primary, clusters));
  }

  return clusters;
}

RawCluster* CaloRawClusterEval::best_cluster_from(PHG4Shower* primary)
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

  if (!get_truth_eval()->is_primary(primary))
  {
    return nullptr;
  }

  primary = get_truth_eval()->get_primary_shower(primary);

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
    std::map<PHG4Shower*, RawCluster*>::iterator iter =
        _cache_best_cluster_from_primary_shower.find(primary);
    if (iter != _cache_best_cluster_from_primary_shower.end())
    {
      return iter->second;
    }
  }

  RawCluster* best_cluster = nullptr;
  float best_energy = FLT_MAX * -1.0;
  std::set<RawCluster*> clusters = all_clusters_from(primary);
  for (auto cluster : clusters)
  {
    if (_strict)
    {
      assert(cluster);
    }
    else if (!cluster)
    {
      ++_errors;
      continue;
    }

    float energy = get_energy_contribution(cluster, primary);
    if (isnan(energy))
    {
      continue;
    }
    if (energy > best_energy)
    {
      best_cluster = cluster;
      best_energy = energy;
    }
  }

  if (_do_cache)
  {
    _cache_best_cluster_from_primary_shower.insert(std::make_pair(primary, best_cluster));
  }

  return best_cluster;
}

float CaloRawClusterEval::get_energy_contribution(RawCluster* cluster, PHG4Shower* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(cluster);
    assert(primary);
  }
  else if (!cluster || !primary)
  {
    ++_errors;
    return NAN;
  }

  if (!get_truth_eval()->is_primary(primary))
  {
    return NAN;
  }

  // reduce cache misses by using only pointer from PrimaryMap
  primary = get_truth_eval()->get_primary_shower(primary);

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
    std::map<std::pair<RawCluster*, PHG4Shower*>, float>::iterator iter =
        _cache_get_energy_contribution_primary_shower.find(std::make_pair(cluster, primary));
    if (iter != _cache_get_energy_contribution_primary_shower.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;

  // loop over all the clustered towers
  if (_usetowerinfo)
  {
    const RawCluster::TowerMap& tower_map = cluster->get_towermap();
    for (auto tower_iter : tower_map)
    {
      RawTowerDefs::keytype tower_key = tower_iter.first;
      unsigned int towerinfo_key = get_towerinfo_key(tower_key);
      TowerInfo* tower = _towerinfos->get_tower_at_key(towerinfo_key);
      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      float edep = get_rawtower_eval()->get_energy_contribution(tower, primary);
      if (!isnan(edep))
      {
        energy += edep;
      }
    }
  }

  else
  {
    RawCluster::TowerConstRange begin_end = cluster->get_towers();
    for (RawCluster::TowerConstIterator iter = begin_end.first;
         iter != begin_end.second;
         ++iter)
    {
      RawTower* tower = _towers->getTower(iter->first);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      float edep = get_rawtower_eval()->get_energy_contribution(tower, primary);
      if (!isnan(edep))
      {
        energy += edep;
      }
    }
  }

  if (_do_cache)
  {
    _cache_get_energy_contribution_primary_shower.insert(std::make_pair(std::make_pair(cluster, primary), energy));
  }

  return energy;
}

std::set<PHG4Particle*> CaloRawClusterEval::all_truth_primary_particles(RawCluster* cluster)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<RawCluster*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_primary_particles.find(cluster);
    if (iter != _cache_all_truth_primary_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_primary_particles;

  std::set<PHG4Shower*> primary_showers = all_truth_primary_showers(cluster);

  for (auto shower : primary_showers)
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

    PHG4Particle* particle = get_truth_eval()->get_primary_particle(shower);

    if (_strict)
    {
      assert(particle);
    }
    else if (!particle)
    {
      ++_errors;
      continue;
    }

    truth_primary_particles.insert(particle);
  }

  if (_do_cache)
  {
    _cache_all_truth_primary_particles.insert(std::make_pair(cluster, truth_primary_particles));
  }

  return truth_primary_particles;
}

PHG4Particle* CaloRawClusterEval::max_truth_primary_particle_by_energy(RawCluster* cluster)
{
  if (!has_reduced_node_pointers())
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
    std::map<RawCluster*, PHG4Particle*>::iterator iter =
        _cache_max_truth_primary_particle_by_energy.find(cluster);
    if (iter != _cache_max_truth_primary_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  PHG4Particle* max_primary = nullptr;
  PHG4Shower* max_shower = max_truth_primary_shower_by_energy(cluster);

  if (max_shower)
  {
    max_primary = get_truth_eval()->get_primary_particle(max_shower);
  }

  if (_do_cache)
  {
    _cache_max_truth_primary_particle_by_energy.insert(std::make_pair(cluster, max_primary));
  }

  return max_primary;
}

std::set<RawCluster*> CaloRawClusterEval::all_clusters_from(PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return std::set<RawCluster*>();
  }

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<RawCluster*>();
  }

  if (!get_truth_eval()->is_primary(primary))
  {
    return std::set<RawCluster*>();
  }

  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++_errors;
    return std::set<RawCluster*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<RawCluster*> >::iterator iter =
        _cache_all_clusters_from_primary_particle.find(primary);
    if (iter != _cache_all_clusters_from_primary_particle.end())
    {
      return iter->second;
    }
  }

  std::set<RawCluster*> clusters;

  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);

  if (shower)
  {
    clusters = all_clusters_from(shower);
  }

  if (_do_cache)
  {
    _cache_all_clusters_from_primary_particle.insert(std::make_pair(primary, clusters));
  }

  return clusters;
}

RawCluster* CaloRawClusterEval::best_cluster_from(PHG4Particle* primary)
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

  if (!get_truth_eval()->is_primary(primary))
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
    std::map<PHG4Particle*, RawCluster*>::iterator iter =
        _cache_best_cluster_from_primary_particle.find(primary);
    if (iter != _cache_best_cluster_from_primary_particle.end())
    {
      return iter->second;
    }
  }

  RawCluster* best_cluster = nullptr;

  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);
  if (shower)
  {
    best_cluster = best_cluster_from(shower);
  }

  if (_do_cache)
  {
    _cache_best_cluster_from_primary_particle.insert(std::make_pair(primary, best_cluster));
  }

  return best_cluster;
}

float CaloRawClusterEval::get_energy_contribution(RawCluster* cluster, PHG4Particle* primary)
{
  if (!has_reduced_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(cluster);
    assert(primary);
  }
  else if (!cluster || !primary)
  {
    ++_errors;
    return NAN;
  }

  if (!get_truth_eval()->is_primary(primary))
  {
    return NAN;
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
    return 0.;
  }

  if (_do_cache)
  {
    std::map<std::pair<RawCluster*, PHG4Particle*>, float>::iterator iter =
        _cache_get_energy_contribution_primary_particle.find(std::make_pair(cluster, primary));
    if (iter != _cache_get_energy_contribution_primary_particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;
  PHG4Shower* shower = get_truth_eval()->get_primary_shower(primary);
  if (shower)
  {
    float edep = get_energy_contribution(cluster, shower);
    if (!isnan(edep))
    {
      energy += edep;
    }
  }

  if (_do_cache)
  {
    _cache_get_energy_contribution_primary_particle.insert(std::make_pair(std::make_pair(cluster, primary), energy));
  }

  return energy;
}

bool CaloRawClusterEval::has_full_node_pointers()
{
  if (!get_rawtower_eval()->has_full_node_pointers())
  {
    return false;
  }

  if (_strict)
  {
    assert(_clusters);
  }
  else if (!_clusters)
  {
    return false;
  }

  if (_strict)
  {
    assert(_towers);
  }
  else if (!_towers)
  {
    return false;
  }

  return true;
}

std::set<PHG4Hit*> CaloRawClusterEval::all_truth_hits(RawCluster* cluster)
{
  if (!has_full_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<RawCluster*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(cluster);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the clustered towers
  if (_usetowerinfo)
  {
    const RawCluster::TowerMap& tower_map = cluster->get_towermap();
    for (auto tower_iter : tower_map)
    {
      RawTowerDefs::keytype tower_key = tower_iter.first;
      unsigned int towerinfo_key = get_towerinfo_key(tower_key);
      TowerInfo* tower = _towerinfos->get_tower_at_key(towerinfo_key);
      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      std::set<PHG4Hit*> new_hits = get_rawtower_eval()->all_truth_hits(tower);

      for (auto g4hit : new_hits)
      {
        if (_strict)
        {
          assert(g4hit);
        }
        else if (!g4hit)
        {
          ++_errors;
          continue;
        }

        truth_hits.insert(g4hit);
      }
    }
  }
  else
  {
    RawCluster::TowerConstRange begin_end = cluster->get_towers();
    for (RawCluster::TowerConstIterator iter = begin_end.first;
         iter != begin_end.second;
         ++iter)
    {
      RawTower* tower = _towers->getTower(iter->first);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      std::set<PHG4Hit*> new_hits = get_rawtower_eval()->all_truth_hits(tower);

      for (auto g4hit : new_hits)
      {
        if (_strict)
        {
          assert(g4hit);
        }
        else if (!g4hit)
        {
          ++_errors;
          continue;
        }

        truth_hits.insert(g4hit);
      }
    }
  }
  if (_do_cache)
  {
    _cache_all_truth_hits.insert(std::make_pair(cluster, truth_hits));
  }

  return truth_hits;
}

void CaloRawClusterEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  std::string nodename = "CLUSTERINFO_" + _caloname;
  if (!_usetowerinfo)
  {
    nodename = "CLUSTER_" + _caloname;
  }
  _clusters = findNode::getClass<RawClusterContainer>(topNode, nodename.c_str());

  std::string towername = "TOWER_CALIB_" + _caloname;
  _towers = findNode::getClass<RawTowerContainer>(topNode, towername.c_str());

  std::string towerinfoname = "TOWERINFO_" + _caloname;
  _towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towerinfoname.c_str());

  return;
}

unsigned int CaloRawClusterEval::get_towerinfo_key(RawTowerDefs::keytype tower_key)
{
  int ix = RawTowerDefs::decode_index2(tower_key);  // iphi
  int iy = RawTowerDefs::decode_index1(tower_key);  // ieta
  RawTowerDefs::CalorimeterId caloid =
      RawTowerDefs::decode_caloid(tower_key);
  // the encoding for calo are actually all the same
  //  this is for safety and furture compatibility(s.l.)
  unsigned int towerinfokey = UINT_MAX;
  if (caloid == RawTowerDefs::CalorimeterId::CEMC)
  {
    towerinfokey = TowerInfoDefs::encode_emcal(iy, ix);
  }
  else if (caloid == RawTowerDefs::CalorimeterId::HCALIN || caloid == RawTowerDefs::CalorimeterId::HCALOUT)
  {
    towerinfokey = TowerInfoDefs::encode_hcal(iy, ix);
  }
  else
  {
    std::cout << "CaloRawClusterEval::get_towerinfo_key - unknown caloid: " << caloid << std::endl;
  }

  return towerinfokey;
}
