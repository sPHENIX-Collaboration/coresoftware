#include "JetRecoEval.h"

#include "CaloEvalStack.h"
#include "CaloRawClusterEval.h"
#include "CaloRawTowerEval.h"  // for CaloRawTowerEval
#include "JetTruthEval.h"
#include "SvtxEvalStack.h"
#include "SvtxTrackEval.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <g4main/PHG4Particle.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>

JetRecoEval::JetRecoEval(PHCompositeNode* topNode,
                         const std::string& recojetname,
                         const std::string& truthjetname)
  : _jettrutheval(topNode, truthjetname)
  , _recojetname(recojetname)
  , _truthjetname(truthjetname)
{
  get_node_pointers(topNode);
}

JetRecoEval::~JetRecoEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      std::cout << "JetRecoEval::~JetRecoEval() - Error Count: " << _errors << std::endl;
    }
  }
}

void JetRecoEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_showers.clear();
  _cache_all_truth_particles.clear();
  _cache_all_truth_jets.clear();
  _cache_max_truth_jet_by_energy.clear();
  _cache_all_jets_from.clear();
  _cache_best_jet_from.clear();
  _cache_get_energy_contribution.clear();
  _cache_get_energy_contribution_src.clear();
  _cache_all_truth_hits.clear();

  _jettrutheval.next_event(topNode);

  get_node_pointers(topNode);
}

void JetRecoEval::set_track_nodename(const std::string& name)
{
  m_TrackNodeName = name;
  _jettrutheval.set_track_nodename(name);
}

std::set<PHG4Shower*> JetRecoEval::all_truth_showers(Jet* recojet)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<PHG4Shower*> >::iterator iter =
        _cache_all_truth_showers.find(recojet);
    if (iter != _cache_all_truth_showers.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Shower*> truth_showers;

  // loop over all the jet constituents, backtrack each reco object to the
  // truth hits and combine with other consituents
  for (auto jter : recojet->get_comp_vec())
  {
      Jet::SRC source = jter.first;
      unsigned int index = jter.second;

    std::set<PHG4Shower*> new_showers;

    if (source == Jet::TRACK)
    {
      if (!_trackmap)
      {
        std::cout << PHWHERE << "ERROR: can't find SvtxTrackMap" << std::endl;
        exit(-1);
      }

      // SvtxTrack* track = _trackmap->get(index);

      // if (_strict) {assert(track);}
      // else if (!track) {++_errors; continue;}

      // new_showers = get_svtx_eval_stack()->get_track_eval()->all_truth_showers(track);
    }

    else if (source == Jet::CEMC_TOWER)
    {
      if (!_cemctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_CEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _cemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_showers = get_cemc_eval_stack()->get_rawtower_eval()->all_truth_primary_showers(tower);
    }
    else if (source == Jet::CEMC_CLUSTER)
    {
      if (!_cemcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_CEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _cemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_showers = get_cemc_eval_stack()->get_rawcluster_eval()->all_truth_primary_showers(cluster);
    }
    else if (source == Jet::EEMC_TOWER)
    {
      if (!_eemctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_EEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _eemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_showers = get_eemc_eval_stack()->get_rawtower_eval()->all_truth_primary_showers(tower);
    }
    else if (source == Jet::EEMC_CLUSTER)
    {
      if (!_eemcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_EEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _eemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_showers = get_eemc_eval_stack()->get_rawcluster_eval()->all_truth_primary_showers(cluster);
    }
    else if (source == Jet::HCALIN_TOWER)
    {
      if (!_hcalintowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_HCALIN" << std::endl;
        exit(-1);
      }

      RawTower* tower = _hcalintowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_showers = get_hcalin_eval_stack()->get_rawtower_eval()->all_truth_primary_showers(tower);
    }
    else if (source == Jet::HCALIN_CLUSTER)
    {
      if (!_hcalinclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_HCALIN" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _hcalinclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_showers = get_hcalin_eval_stack()->get_rawcluster_eval()->all_truth_primary_showers(cluster);
    }
    else if (source == Jet::HCALOUT_TOWER)
    {
      if (!_hcalouttowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_HCALOUT" << std::endl;
        exit(-1);
      }

      RawTower* tower = _hcalouttowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_showers = get_hcalout_eval_stack()->get_rawtower_eval()->all_truth_primary_showers(tower);
    }
    else if (source == Jet::HCALOUT_CLUSTER)
    {
      if (!_hcaloutclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_HCALOUT" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _hcaloutclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_showers = get_hcalout_eval_stack()->get_rawcluster_eval()->all_truth_primary_showers(cluster);
    }
    else if (source == Jet::FEMC_TOWER)
    {
      if (!_femctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_FEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _femctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_showers = get_femc_eval_stack()->get_rawtower_eval()->all_truth_primary_showers(tower);
    }
    else if (source == Jet::FEMC_CLUSTER)
    {
      if (!_femcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_FEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _femcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_showers = get_femc_eval_stack()->get_rawcluster_eval()->all_truth_primary_showers(cluster);
    }
    else if (source == Jet::FHCAL_TOWER)
    {
      if (!_fhcaltowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_FHCAL" << std::endl;
        exit(-1);
      }

      RawTower* tower = _fhcaltowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_showers = get_fhcal_eval_stack()->get_rawtower_eval()->all_truth_primary_showers(tower);
    }
    else if (source == Jet::FHCAL_CLUSTER)
    {
      if (!_fhcalclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_FHCAL" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _fhcalclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_showers = get_fhcal_eval_stack()->get_rawcluster_eval()->all_truth_primary_showers(cluster);
    }

    for (auto new_shower : new_showers)
    {
      truth_showers.insert(new_shower);
    }
  }

  if (_do_cache)
  {
    _cache_all_truth_showers.insert(std::make_pair(recojet, truth_showers));
  }

  return truth_showers;
}

std::set<PHG4Particle*> JetRecoEval::all_truth_particles(Jet* recojet)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(recojet);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  // loop over all the jet constituents, backtrack each reco object to the
  // truth hits and combine with other consituents
  for (auto jter : recojet->get_comp_vec())
  {
      Jet::SRC source = jter.first;
      unsigned int index = jter.second;

    std::set<PHG4Particle*> new_particles;

    if (source == Jet::TRACK)
    {
      if (!_trackmap)
      {
        std::cout << PHWHERE << "ERROR: can't find TrackMap" << std::endl;
        exit(-1);
      }

      SvtxTrack* track = _trackmap->get(index);

      if (_strict)
      {
        assert(track);
      }
      else if (!track)
      {
        ++_errors;
        continue;
      }

      new_particles = get_svtx_eval_stack()->get_track_eval()->all_truth_particles(track);
    }
    else if (source == Jet::CEMC_TOWER)
    {
      if (!_cemctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_CEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _cemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_particles = get_cemc_eval_stack()->get_rawtower_eval()->all_truth_primary_particles(tower);
    }
    else if (source == Jet::CEMC_CLUSTER)
    {
      if (!_cemcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_CEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _cemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_particles = get_cemc_eval_stack()->get_rawcluster_eval()->all_truth_primary_particles(cluster);
    }
    else if (source == Jet::EEMC_TOWER)
    {
      if (!_eemctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_EEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _eemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_particles = get_eemc_eval_stack()->get_rawtower_eval()->all_truth_primary_particles(tower);
    }
    else if (source == Jet::EEMC_CLUSTER)
    {
      if (!_eemcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_EEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _eemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_particles = get_eemc_eval_stack()->get_rawcluster_eval()->all_truth_primary_particles(cluster);
    }
    else if (source == Jet::HCALIN_TOWER)
    {
      if (!_hcalintowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_HCALIN" << std::endl;
        exit(-1);
      }

      RawTower* tower = _hcalintowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_particles = get_hcalin_eval_stack()->get_rawtower_eval()->all_truth_primary_particles(tower);
    }
    else if (source == Jet::HCALIN_CLUSTER)
    {
      if (!_hcalinclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_HCALIN" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _hcalinclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_particles = get_hcalin_eval_stack()->get_rawcluster_eval()->all_truth_primary_particles(cluster);
    }
    else if (source == Jet::HCALOUT_TOWER)
    {
      if (!_hcalouttowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_HCALOUT" << std::endl;
        exit(-1);
      }

      RawTower* tower = _hcalouttowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_particles = get_hcalout_eval_stack()->get_rawtower_eval()->all_truth_primary_particles(tower);
    }
    else if (source == Jet::HCALOUT_CLUSTER)
    {
      if (!_hcaloutclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_HCALOUT" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _hcaloutclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_particles = get_hcalout_eval_stack()->get_rawcluster_eval()->all_truth_primary_particles(cluster);
    }
    else if (source == Jet::FEMC_TOWER)
    {
      if (!_femctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_FEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _femctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_particles = get_femc_eval_stack()->get_rawtower_eval()->all_truth_primary_particles(tower);
    }
    else if (source == Jet::FEMC_CLUSTER)
    {
      if (!_femcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_FEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _femcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_particles = get_femc_eval_stack()->get_rawcluster_eval()->all_truth_primary_particles(cluster);
    }
    else if (source == Jet::FHCAL_TOWER)
    {
      if (!_fhcaltowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_FHCAL" << std::endl;
        exit(-1);
      }

      RawTower* tower = _fhcaltowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_particles = get_fhcal_eval_stack()->get_rawtower_eval()->all_truth_primary_particles(tower);
    }
    else if (source == Jet::FHCAL_CLUSTER)
    {
      if (!_fhcalclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_FHCAL" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _fhcalclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_particles = get_fhcal_eval_stack()->get_rawcluster_eval()->all_truth_primary_particles(cluster);
    }

    for (auto new_particle : new_particles)
    {
      truth_particles.insert(new_particle);
    }
  }

  if (_do_cache)
  {
    _cache_all_truth_particles.insert(std::make_pair(recojet, truth_particles));
  }

  return truth_particles;
}

std::set<Jet*> JetRecoEval::all_truth_jets(Jet* recojet)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return std::set<Jet*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<Jet*> >::iterator iter =
        _cache_all_truth_jets.find(recojet);
    if (iter != _cache_all_truth_jets.end())
    {
      return iter->second;
    }
  }

  std::set<Jet*> truth_jets;

  // get all truth particles (this can include muons and other truth excludes)...
  std::set<PHG4Particle*> particles = all_truth_particles(recojet);

  // backtrack from the truth particles to the truth jets...
  for (auto particle : particles)
  {
    if (_strict)
    {
      assert(particle);
    }
    else if (!particle)
    {
      ++_errors;
      continue;
    }

    Jet* truth_jet = _jettrutheval.get_truth_jet(particle);
    if (!truth_jet)
    {
      continue;
    }

    truth_jets.insert(truth_jet);
  }

  if (_do_cache)
  {
    _cache_all_truth_jets.insert(std::make_pair(recojet, truth_jets));
  }

  return truth_jets;
}

Jet* JetRecoEval::max_truth_jet_by_energy(Jet* recojet)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<Jet*, Jet*>::iterator iter =
        _cache_max_truth_jet_by_energy.find(recojet);
    if (iter != _cache_max_truth_jet_by_energy.end())
    {
      return iter->second;
    }
  }

  Jet* truthjet = nullptr;
  float max_energy = FLT_MAX * -1.0;

  std::set<Jet*> truthjets = all_truth_jets(recojet);
  for (auto candidate : truthjets)
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

    float energy = get_energy_contribution(recojet, candidate);
    if (energy > max_energy)
    {
      truthjet = candidate;
      max_energy = energy;
    }
  }

  if (_do_cache)
  {
    _cache_max_truth_jet_by_energy.insert(std::make_pair(recojet, truthjet));
  }

  return truthjet;
}

std::set<Jet*> JetRecoEval::all_jets_from(Jet* truthjet)
{
  if (_strict)
  {
    assert(truthjet);
  }
  else if (!truthjet)
  {
    ++_errors;
    return std::set<Jet*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<Jet*> >::iterator iter =
        _cache_all_jets_from.find(truthjet);
    if (iter != _cache_all_jets_from.end())
    {
      return iter->second;
    }
  }

  std::set<Jet*> recojets;

  // loop over all reco jets
  for (auto recojet : *_recojets)
  {
    /* Jet* recojet = _recojet.second; */

    // if this jet back tracks to the truth jet
    std::set<Jet*> truthcandidates = all_truth_jets(recojet);
    for (auto truthcandidate : truthcandidates)
    {
      if (_strict)
      {
        assert(truthcandidate);
      }
      else if (!truthcandidate)
      {
        ++_errors;
        continue;
      }

      if (truthcandidate->get_id() == truthjet->get_id())
      {
        recojets.insert(recojet);
      }
    }
  }

  if (_do_cache)
  {
    _cache_all_jets_from.insert(std::make_pair(truthjet, recojets));
  }

  return recojets;
}

Jet* JetRecoEval::best_jet_from(Jet* truthjet)
{

  if (_strict)
  {
    assert(truthjet);
  }
  else if (!truthjet)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<Jet*, Jet*>::iterator iter =
        _cache_best_jet_from.find(truthjet);
    if (iter != _cache_best_jet_from.end())
    {
      return iter->second;
    }
  }

  Jet* bestrecojet = nullptr;
  float max_energy = FLT_MAX * -1.0;

  std::set<Jet*> recojets = all_jets_from(truthjet);
  for (auto recojet : recojets)
  {
    if (_strict)
    {
      assert(recojet);
    }
    else if (!recojet)
    {
      ++_errors;
      continue;
    }

    float energy = get_energy_contribution(recojet, truthjet);
    if (energy > max_energy)
    {
      bestrecojet = recojet;
      max_energy = energy;
    }
  }

  if (_do_cache)
  {
    _cache_best_jet_from.insert(std::make_pair(truthjet, bestrecojet));
  }
  return bestrecojet;
}

Jet* JetRecoEval::unique_reco_jet_from_truth(Jet* truthjet)
{
  if (_strict)
  {
    assert(truthjet);
  }
  else if (!truthjet)
  {
    ++_errors;
    return nullptr;
  }

  Jet* recojet = best_jet_from(truthjet);

  if (recojet)
  {
    Jet* back_matching = max_truth_jet_by_energy(recojet);

    if (back_matching->get_id() == truthjet->get_id())
    {
      return recojet;  // uniquely matched
    }
    else
    {
      return nullptr;
    }
  }
  else
  {
    return nullptr;
  }
}

Jet* JetRecoEval::unique_truth_jet_from_reco(Jet* recojet)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return nullptr;
  }

  Jet* truthjet = max_truth_jet_by_energy(recojet);

  if (truthjet)
  {
    Jet* back_matching = best_jet_from(truthjet);

    if (back_matching->get_id() == recojet->get_id())
    {
      return truthjet;  // uniquely matched
    }
    else
    {
      return nullptr;
    }
  }
  else
  {
    return nullptr;
  }
}

// overlap calculations
float JetRecoEval::get_energy_contribution(Jet* recojet, Jet* truthjet)
{
  if (_strict)
  {
    assert(recojet);
    assert(truthjet);
  }
  else if (!recojet || !truthjet)
  {
    ++_errors;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<std::pair<Jet*, Jet*>, float>::iterator iter =
        _cache_get_energy_contribution.find(std::make_pair(recojet, truthjet));
    if (iter != _cache_get_energy_contribution.end())
    {
      return iter->second;
    }
  }

  float energy_contribution = 0.0;

  std::set<PHG4Particle*> truthjetcomp = get_truth_eval()->all_truth_particles(truthjet);
  // loop over all truthjet constituents
  for (auto truthparticle : truthjetcomp)
  {
    if (_strict)
    {
      assert(truthparticle);
    }
    else if (!truthparticle)
    {
      ++_errors;
      continue;
    }

    for (auto jter : recojet->get_comp_vec())
    {
      Jet::SRC source = jter.first;
      unsigned int index = jter.second;

      float energy = 0.0;

      if (source == Jet::TRACK)
      {
        SvtxTrack* track = _trackmap->get(index);

        if (_strict)
        {
          assert(track);
        }
        else if (!track)
        {
          ++_errors;
          continue;
        }

        PHG4Particle* maxtruthparticle = get_svtx_eval_stack()->get_track_eval()->max_truth_particle_by_nclusters(track);

        if (maxtruthparticle == nullptr)
        {
          // in extreme rare cases, noise hits can make a track with no maxtruthparticle matched
          energy = 0;
        }
        else if (maxtruthparticle->get_track_id() == truthparticle->get_track_id())
        {
          energy = track->get_p();
        }
      }
      else if (source == Jet::CEMC_TOWER)
      {
        RawTower* tower = _cemctowers->getTower(index);

        if (_strict)
        {
          assert(tower);
        }
        else if (!tower)
        {
          ++_errors;
          continue;
        }

        energy = get_cemc_eval_stack()->get_rawtower_eval()->get_energy_contribution(tower, truthparticle);
      }
      else if (source == Jet::CEMC_CLUSTER)
      {
        RawCluster* cluster = _cemcclusters->getCluster(index);

        if (_strict)
        {
          assert(cluster);
        }
        else if (!cluster)
        {
          ++_errors;
          continue;
        }

        energy = get_cemc_eval_stack()->get_rawcluster_eval()->get_energy_contribution(cluster, truthparticle);
      }
      else if (source == Jet::EEMC_TOWER)
      {
        RawTower* tower = _eemctowers->getTower(index);

        if (_strict)
        {
          assert(tower);
        }
        else if (!tower)
        {
          ++_errors;
          continue;
        }

        energy = get_eemc_eval_stack()->get_rawtower_eval()->get_energy_contribution(tower, truthparticle);
      }
      else if (source == Jet::EEMC_CLUSTER)
      {
        RawCluster* cluster = _eemcclusters->getCluster(index);

        if (_strict)
        {
          assert(cluster);
        }
        else if (!cluster)
        {
          ++_errors;
          continue;
        }

        energy = get_eemc_eval_stack()->get_rawcluster_eval()->get_energy_contribution(cluster, truthparticle);
      }
      else if (source == Jet::HCALIN_TOWER)
      {
        RawTower* tower = _hcalintowers->getTower(index);

        if (_strict)
        {
          assert(tower);
        }
        else if (!tower)
        {
          ++_errors;
          continue;
        }

        energy = get_hcalin_eval_stack()->get_rawtower_eval()->get_energy_contribution(tower, truthparticle);
      }
      else if (source == Jet::HCALIN_CLUSTER)
      {
        RawCluster* cluster = _hcalinclusters->getCluster(index);

        if (_strict)
        {
          assert(cluster);
        }
        else if (!cluster)
        {
          ++_errors;
          continue;
        }

        energy = get_hcalin_eval_stack()->get_rawcluster_eval()->get_energy_contribution(cluster, truthparticle);
      }
      else if (source == Jet::HCALOUT_TOWER)
      {
        RawTower* tower = _hcalouttowers->getTower(index);

        if (_strict)
        {
          assert(tower);
        }
        else if (!tower)
        {
          ++_errors;
          continue;
        }

        energy = get_hcalout_eval_stack()->get_rawtower_eval()->get_energy_contribution(tower, truthparticle);
      }
      else if (source == Jet::HCALOUT_CLUSTER)
      {
        RawCluster* cluster = _hcaloutclusters->getCluster(index);

        if (_strict)
        {
          assert(cluster);
        }
        else if (!cluster)
        {
          ++_errors;
          continue;
        }

        energy = get_hcalout_eval_stack()->get_rawcluster_eval()->get_energy_contribution(cluster, truthparticle);
      }
      else if (source == Jet::FEMC_TOWER)
      {
        RawTower* tower = _femctowers->getTower(index);

        if (_strict)
        {
          assert(tower);
        }
        else if (!tower)
        {
          ++_errors;
          continue;
        }

        energy = get_femc_eval_stack()->get_rawtower_eval()->get_energy_contribution(tower, truthparticle);
      }
      else if (source == Jet::FEMC_CLUSTER)
      {
        RawCluster* cluster = _femcclusters->getCluster(index);

        if (_strict)
        {
          assert(cluster);
        }
        else if (!cluster)
        {
          ++_errors;
          continue;
        }

        energy = get_femc_eval_stack()->get_rawcluster_eval()->get_energy_contribution(cluster, truthparticle);
      }
      else if (source == Jet::FHCAL_TOWER)
      {
        RawTower* tower = _fhcaltowers->getTower(index);

        if (_strict)
        {
          assert(tower);
        }
        else if (!tower)
        {
          ++_errors;
          continue;
        }

        energy = get_fhcal_eval_stack()->get_rawtower_eval()->get_energy_contribution(tower, truthparticle);
      }
      else if (source == Jet::FHCAL_CLUSTER)
      {
        RawCluster* cluster = _fhcalclusters->getCluster(index);

        if (_strict)
        {
          assert(cluster);
        }
        else if (!cluster)
        {
          ++_errors;
          continue;
        }

        energy = get_fhcal_eval_stack()->get_rawcluster_eval()->get_energy_contribution(cluster, truthparticle);
      }

      energy_contribution += energy;
    }
  }

  if (_do_cache)
  {
    _cache_get_energy_contribution.insert(std::make_pair(std::make_pair(recojet, truthjet), energy_contribution));
  }

  return energy_contribution;
}

// overlap calculations
float JetRecoEval::get_energy_contribution(Jet* recojet, Jet::SRC src)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<std::pair<Jet*, Jet::SRC>, float>::iterator iter =
        _cache_get_energy_contribution_src.find(std::make_pair(recojet, src));
    if (iter != _cache_get_energy_contribution_src.end())
    {
      return iter->second;
    }
  }
  
  float energy = 0.0;

  // loop over all recojet constituents
  for (Jet::ITER_comp_vec jter = recojet->comp_begin(src);
       jter != recojet->comp_end(src);
       ++jter)
  {
    Jet::SRC source = jter->first;
    assert(source == src);  // jet container consistency check
    unsigned int index = jter->second;

    if (source == Jet::TRACK)
    {
      SvtxTrack* track = _trackmap->get(index);
      energy += track->get_p();
    }
    else if (source == Jet::CEMC_TOWER)
    {
      RawTower* tower = _cemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      energy += tower->get_energy();
    }
    else if (source == Jet::CEMC_CLUSTER)
    {
      RawCluster* cluster = _cemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      energy += cluster->get_energy();
    }
    else if (source == Jet::EEMC_TOWER)
    {
      RawTower* tower = _eemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      energy += tower->get_energy();
    }
    else if (source == Jet::EEMC_CLUSTER)
    {
      RawCluster* cluster = _eemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      energy += cluster->get_energy();
    }
    else if (source == Jet::HCALIN_TOWER)
    {
      RawTower* tower = _hcalintowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      energy += tower->get_energy();
    }
    else if (source == Jet::HCALIN_CLUSTER)
    {
      RawCluster* cluster = _hcalinclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      energy += cluster->get_energy();
    }
    else if (source == Jet::HCALOUT_TOWER)
    {
      RawTower* tower = _hcalouttowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      energy += tower->get_energy();
    }
    else if (source == Jet::HCALOUT_CLUSTER)
    {
      RawCluster* cluster = _hcaloutclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      energy += cluster->get_energy();
    }
    else if (source == Jet::FEMC_TOWER)
    {
      RawTower* tower = _femctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      energy += tower->get_energy();
    }
    else if (source == Jet::FEMC_CLUSTER)
    {
      RawCluster* cluster = _femcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      energy += cluster->get_energy();
    }
    else if (source == Jet::FHCAL_TOWER)
    {
      RawTower* tower = _fhcaltowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      energy += tower->get_energy();
    }
    else if (source == Jet::FHCAL_CLUSTER)
    {
      RawCluster* cluster = _fhcalclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      energy += cluster->get_energy();
    }  // else if (source == Jet::FHCAL_CLUSTER)

  }  // for (Jet::ConstIter jter = recojet->lower_bound_comp(src);

  if (_do_cache)
  {
    _cache_get_energy_contribution_src.insert(std::make_pair(std::make_pair(recojet, src), energy));
  }

  return energy;
}

std::set<PHG4Hit*> JetRecoEval::all_truth_hits(Jet* recojet)
{
  if (_strict)
  {
    assert(recojet);
  }
  else if (!recojet)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(recojet);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the jet constituents, backtrack each reco object to the
  // truth hits and combine with other consituents

  for (auto jter : recojet->get_comp_vec())
  {
      Jet::SRC source = jter.first;
      unsigned int index = jter.second;

    std::set<PHG4Hit*> new_hits;

    if (source == Jet::TRACK)
    {
      if (!_trackmap)
      {
        std::cout << PHWHERE << "ERROR: can't find SvtxTrackMap" << std::endl;
        exit(-1);
      }

      SvtxTrack* track = _trackmap->get(index);

      if (_strict)
      {
        assert(track);
      }
      else if (!track)
      {
        ++_errors;
        continue;
      }

      new_hits = get_svtx_eval_stack()->get_track_eval()->all_truth_hits(track);
    }
    else if (source == Jet::CEMC_TOWER)
    {
      if (!_cemctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_CEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _cemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_hits = get_cemc_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);
    }
    else if (source == Jet::CEMC_CLUSTER)
    {
      if (!_cemcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_CEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _cemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_hits = get_cemc_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);
    }
    else if (source == Jet::EEMC_TOWER)
    {
      if (!_eemctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_EEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _eemctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_hits = get_eemc_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);
    }
    else if (source == Jet::EEMC_CLUSTER)
    {
      if (!_eemcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_EEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _eemcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_hits = get_eemc_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);
    }
    else if (source == Jet::HCALIN_TOWER)
    {
      if (!_hcalintowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_HCALIN" << std::endl;
        exit(-1);
      }

      RawTower* tower = _hcalintowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_hits = get_hcalin_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);
    }
    else if (source == Jet::HCALIN_CLUSTER)
    {
      if (!_hcalinclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_HCALIN" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _hcalinclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_hits = get_hcalin_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);
    }
    else if (source == Jet::HCALOUT_TOWER)
    {
      if (!_hcalouttowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_HCALOUT" << std::endl;
        exit(-1);
      }

      RawTower* tower = _hcalouttowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_hits = get_hcalout_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);
    }
    else if (source == Jet::HCALOUT_CLUSTER)
    {
      if (!_hcaloutclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_HCALOUT" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _hcaloutclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_hits = get_hcalout_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);
    }
    else if (source == Jet::FEMC_TOWER)
    {
      if (!_femctowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_FEMC" << std::endl;
        exit(-1);
      }

      RawTower* tower = _femctowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_hits = get_femc_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);
    }
    else if (source == Jet::FEMC_CLUSTER)
    {
      if (!_femcclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_FEMC" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _femcclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_hits = get_femc_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);
    }
    else if (source == Jet::FHCAL_TOWER)
    {
      if (!_fhcaltowers)
      {
        std::cout << PHWHERE << "ERROR: can't find TOWER_FHCAL" << std::endl;
        exit(-1);
      }

      RawTower* tower = _fhcaltowers->getTower(index);

      if (_strict)
      {
        assert(tower);
      }
      else if (!tower)
      {
        ++_errors;
        continue;
      }

      new_hits = get_fhcal_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);
    }
    else if (source == Jet::FHCAL_CLUSTER)
    {
      if (!_fhcalclusters)
      {
        std::cout << PHWHERE << "ERROR: can't find CLUSTER_FHCAL" << std::endl;
        exit(-1);
      }

      RawCluster* cluster = _fhcalclusters->getCluster(index);

      if (_strict)
      {
        assert(cluster);
      }
      else if (!cluster)
      {
        ++_errors;
        continue;
      }

      new_hits = get_fhcal_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);
    }

    for (auto new_hit : new_hits)
    {
      truth_hits.insert(new_hit);
    }
  }

  if (_do_cache)
  {
    _cache_all_truth_hits.insert(std::make_pair(recojet, truth_hits));
  }

  return truth_hits;
}

void JetRecoEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  _recojets = findNode::getClass<JetContainer>(topNode, _recojetname.c_str());
  if (!_recojets)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << _recojetname << std::endl;
    exit(-1);
  }

  _truthjets = findNode::getClass<JetContainer>(topNode, _truthjetname.c_str());
  if (!_truthjets)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << _truthjetname << std::endl;
    exit(-1);
  }

  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_TrackNodeName);
  if (!_trackmap)
  {
    _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
  }
  _cemctowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  _hcalintowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  _hcalouttowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
  _femctowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_FEMC");
  _fhcaltowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_FHCAL");
  _eemctowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_EEMC");
  _cemcclusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
  _hcalinclusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
  _hcaloutclusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
  _femcclusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_FEMC");
  _fhcalclusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_FHCAL");
  _eemcclusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_EEMC");

  return;
}
