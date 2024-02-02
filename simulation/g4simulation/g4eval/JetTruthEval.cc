#include "JetTruthEval.h"

#include "CaloTruthEval.h"
#include "SvtxTruthEval.h"

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

JetTruthEval::JetTruthEval(PHCompositeNode* topNode,
                           const std::string& truthjetname)
  : _truthjetname(truthjetname)
  , _svtxevalstack(topNode)
  , _cemcevalstack(topNode, "CEMC")
  , _hcalinevalstack(topNode, "HCALIN")
  , _hcaloutevalstack(topNode, "HCALOUT")
  , _femcevalstack(topNode, "FEMC")
  , _fhcalevalstack(topNode, "FHCAL")
  , _eemcevalstack(topNode, "EEMC")
{
  get_node_pointers(topNode);
}

JetTruthEval::~JetTruthEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      std::cout << "JetTruthEval::~JetTruthEval() - Error Count: " << _errors << std::endl;
    }
  }
}

void JetTruthEval::next_event(PHCompositeNode* topNode)
{
  _svtxevalstack.next_event(topNode);
  _cemcevalstack.next_event(topNode);
  _hcalinevalstack.next_event(topNode);
  _hcaloutevalstack.next_event(topNode);
  _femcevalstack.next_event(topNode);
  _fhcalevalstack.next_event(topNode);

  _cache_all_truth_particles.clear();
  _cache_all_truth_showers.clear();
  _cache_all_truth_hits.clear();
  _cache_get_truth_jet.clear();

  get_node_pointers(topNode);
}

std::set<PHG4Particle*> JetTruthEval::all_truth_particles(Jet* truthjet)
{
  if (_strict)
  {
    assert(truthjet);
  }
  else if (!truthjet)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(truthjet);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  // loop over all the entries in the truthjet
  /* for (Jet::ConstIter iter = truthjet->comp_begin(); */
  /*      iter != truthjet->comp_end(); */
  /*      ++iter) */
  for (const auto& iter : truthjet->get_comp_vec())  // vector of pair<Jet::SRC, unsigned int>
  {
    Jet::SRC source = iter.first;
    unsigned int index = iter.second;
    if (source != Jet::PARTICLE)
    {
      std::cout << PHWHERE << " truth jet contains something other than particles!" << std::endl;
      exit(-1);
    }

    PHG4Particle* truth_particle = _truthinfo->GetParticle(index);

    if (_strict)
    {
      assert(truth_particle);
    }
    else if (!truth_particle)
    {
      ++_errors;
      continue;
    }

    truth_particles.insert(truth_particle);
  }

  if (_do_cache)
  {
    _cache_all_truth_particles.insert(std::make_pair(truthjet, truth_particles));
  }

  return truth_particles;
}

std::set<PHG4Shower*> JetTruthEval::all_truth_showers(Jet* truthjet)
{
  if (_strict)
  {
    assert(truthjet);
  }
  else if (!truthjet)
  {
    ++_errors;
    return std::set<PHG4Shower*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<PHG4Shower*> >::iterator iter =
        _cache_all_truth_showers.find(truthjet);
    if (iter != _cache_all_truth_showers.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Shower*> truth_showers;

  std::set<PHG4Particle*> truth_particles = all_truth_particles(truthjet);

  // loop over all the entries in the truthjet
  for (auto particle : truth_particles)
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

    // any calo truth eval module would work here...
    CaloTruthEval* cemc_truth_eval = _cemcevalstack.get_truth_eval();
    PHG4Shower* shower = cemc_truth_eval->get_primary_shower(particle);
    if (shower)
    {
      truth_showers.insert(shower);
    }
  }

  if (_do_cache)
  {
    _cache_all_truth_showers.insert(std::make_pair(truthjet, truth_showers));
  }

  return truth_showers;
}

std::set<PHG4Hit*> JetTruthEval::all_truth_hits(Jet* truthjet)
{
  if (_strict)
  {
    assert(truthjet);
  }
  else if (!truthjet)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<Jet*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(truthjet);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  std::set<PHG4Particle*> truth_particles = all_truth_particles(truthjet);

  // loop over all the entries in the truthjet
  for (auto particle : truth_particles)
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

    // ask the svtx truth eval to backtrack the particles to g4hits
    SvtxTruthEval* svtx_truth_eval = _svtxevalstack.get_truth_eval();
    std::set<PHG4Hit*> svtx_g4hits = svtx_truth_eval->all_truth_hits(particle);

    for (auto g4hit : svtx_g4hits)
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

  if (_do_cache)
  {
    _cache_all_truth_hits.insert(std::make_pair(truthjet, truth_hits));
  }

  return truth_hits;
}

Jet* JetTruthEval::get_truth_jet(PHG4Particle* particle)
{
  if (_strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, Jet*>::iterator iter =
        _cache_get_truth_jet.find(particle);
    if (iter != _cache_get_truth_jet.end())
    {
      return iter->second;
    }
  }

  Jet* truth_jet = nullptr;

  // loop over all jets and look for this particle...
  for (auto candidate : *_truthjets)
  {
    /* Jet* candidate = _truthjet.second; */

    // loop over all consituents and look for this particle
    for (const std::pair<Jet::SRC, unsigned int>& jter 
        : candidate->get_comp_vec())
    {
      unsigned int index = jter.second;

      PHG4Particle* constituent = _truthinfo->GetParticle(index);
      if (_strict)
      {
        assert(constituent);
      }
      else if (!constituent)
      {
        ++_errors;
        continue;
      }

      if (get_svtx_eval_stack()->get_truth_eval()->are_same_particle(constituent, particle))
      {
        truth_jet = candidate;
        break;
      }
    }

    if (truth_jet)
    {
      break;
    }
  }

  if (_do_cache)
  {
    _cache_get_truth_jet.insert(std::make_pair(particle, truth_jet));
  }

  return truth_jet;
}

void JetTruthEval::get_node_pointers(PHCompositeNode* topNode)
{
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truthinfo)
  {
    std::cout << PHWHERE << " ERROR: Can't find G4TruthInfo" << std::endl;
    exit(-1);
  }

  _truthjets = findNode::getClass<JetContainer>(topNode, _truthjetname);
  if (!_truthjets)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << _truthjetname << std::endl;
    exit(-1);
  }

  return;
}
