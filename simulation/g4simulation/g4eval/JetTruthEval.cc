#include "JetTruthEval.h"

#include "CaloTruthEval.h"
#include "SvtxTruthEval.h"

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

using namespace std;

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
  , _truthinfo(nullptr)
  , _truthjets(nullptr)
  , _strict(false)
  , _verbosity(1)
  , _errors(0)
  , _do_cache(true)
  , _cache_all_truth_particles()
  , _cache_all_truth_showers()
  , _cache_all_truth_hits()
  , _cache_get_truth_jet()
{
  get_node_pointers(topNode);
}

JetTruthEval::~JetTruthEval()
{
  if (_verbosity > 0)
  {
    if ((_errors > 0) || (_verbosity > 1))
    {
      cout << "JetTruthEval::~JetTruthEval() - Error Count: " << _errors << endl;
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
  for (Jet::ConstIter iter = truthjet->begin_comp();
       iter != truthjet->end_comp();
       ++iter)
  {
    Jet::SRC source = iter->first;
    unsigned int index = iter->second;
    if (source != Jet::PARTICLE)
    {
      cout << PHWHERE << " truth jet contains something other than particles!" << endl;
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

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(truthjet, truth_particles));

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
  for (std::set<PHG4Particle*>::iterator iter = truth_particles.begin();
       iter != truth_particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;

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

  if (_do_cache) _cache_all_truth_showers.insert(make_pair(truthjet, truth_showers));

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
  for (std::set<PHG4Particle*>::iterator iter = truth_particles.begin();
       iter != truth_particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;

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

    for (std::set<PHG4Hit*>::iterator jter = svtx_g4hits.begin();
         jter != svtx_g4hits.end();
         ++jter)
    {
      PHG4Hit* g4hit = *jter;

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

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(truthjet, truth_hits));

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
  for (JetMap::Iter iter = _truthjets->begin();
       iter != _truthjets->end();
       ++iter)
  {
    Jet* candidate = iter->second;

    // loop over all consituents and look for this particle
    for (Jet::ConstIter jter = candidate->begin_comp();
         jter != candidate->end_comp();
         ++jter)
    {
      unsigned int index = jter->second;

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

    if (truth_jet) break;
  }

  if (_do_cache) _cache_get_truth_jet.insert(make_pair(particle, truth_jet));

  return truth_jet;
}

void JetTruthEval::get_node_pointers(PHCompositeNode* topNode)
{
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truthinfo)
  {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  _truthjets = findNode::getClass<JetMap>(topNode, _truthjetname.c_str());
  if (!_truthjets)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _truthjetname << endl;
    exit(-1);
  }

  return;
}
