#include "BaseTruthEval.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>

using namespace std;

BaseTruthEval::BaseTruthEval(PHCompositeNode* topNode)
  : m_TruthInfo(nullptr)
  , m_Strict(false)
  , m_Verbosity(0)
  , m_Errors(0)
{
  get_node_pointers(topNode);
}

BaseTruthEval::~BaseTruthEval()
{
  if (m_Verbosity > 0)
  {
    if ((m_Errors > 0) || (m_Verbosity > 1))
    {
      cout << "BaseTruthEval::~BaseTruthEval() - Error Count: " << m_Errors << endl;
    }
  }
}

void BaseTruthEval::next_event(PHCompositeNode* topNode)
{
  get_node_pointers(topNode);
}

int BaseTruthEval::get_embed(PHG4Particle* particle)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return 0;
  }

  if (m_Strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++m_Errors;
    return 0;
  }

  //  if (!is_primary(particle)) return 0;
  // allow the secondary particles to be tagged with its embedding ID of tis primary particle

  PHG4Particle* primary = get_primary_particle(particle);
  if (m_Strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++m_Errors;
    return 0;
  }

  return m_TruthInfo->isEmbeded(primary->get_track_id());
}

PHG4VtxPoint* BaseTruthEval::get_vertex(PHG4Particle* particle)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++m_Errors;
    return nullptr;
  }

  PHG4VtxPoint* vtx = m_TruthInfo->GetVtx(particle->get_vtx_id());
  if (m_Strict)
  {
    assert(vtx);
  }
  else if (!vtx)
  {
    m_Errors++;
  }

  return vtx;
}

bool BaseTruthEval::is_primary(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++m_Errors;
    return false;
  }

  bool is_primary = false;
  if (shower->get_parent_shower_id() == 0)
  {
    is_primary = true;
  }

  return is_primary;
}

bool BaseTruthEval::is_primary(PHG4Particle* particle)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++m_Errors;
    return false;
  }

  bool is_primary = false;
  if (particle->get_parent_id() == 0)
  {
    is_primary = true;
  }

  return is_primary;
}

PHG4Shower* BaseTruthEval::get_primary_shower(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++m_Errors;
    return nullptr;
  }

  if (is_primary(shower)) return shower;

  while (!is_primary(shower))
  {
    shower = m_TruthInfo->GetShower(shower->get_parent_shower_id());

    if (m_Strict)
    {
      assert(shower);
    }
    else if (!shower)
    {
      ++m_Errors;
      break;
    }
  }

  return shower;
}

PHG4Shower* BaseTruthEval::get_primary_shower(PHG4Particle* particle)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++m_Errors;
    return nullptr;
  }

  if (!is_primary(particle)) particle = get_primary_particle(particle);

  PHG4Shower* returnval = nullptr;

  PHG4TruthInfoContainer::ShowerRange range = m_TruthInfo->GetPrimaryShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    PHG4Shower* shower = iter->second;
    if (shower->get_parent_particle_id() == particle->get_track_id())
    {
      returnval = shower;
      break;
    }
  }

  return returnval;
}

PHG4Particle* BaseTruthEval::get_primary_particle(PHG4Particle* particle)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++m_Errors;
    return nullptr;
  }

  if (is_primary(particle)) return particle;

  PHG4Particle* returnval = m_TruthInfo->GetPrimaryParticle(particle->get_primary_id());

  if (m_Strict)
  {
    assert(returnval);
  }
  else if (!returnval)
  {
    ++m_Errors;
  }

  return returnval;
}

PHG4Particle* BaseTruthEval::get_primary_particle(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++m_Errors;
    return nullptr;
  }

  PHG4Particle* parent_particle = m_TruthInfo->GetParticle(shower->get_parent_particle_id());

  if (m_Strict)
  {
    assert(parent_particle);
  }
  else if (!parent_particle)
  {
    ++m_Errors;
    return nullptr;
  }

  PHG4Particle* primary_particle = get_primary_particle(parent_particle);

  if (m_Strict)
  {
    assert(primary_particle);
  }
  else if (!primary_particle)
  {
    ++m_Errors;
    return nullptr;
  }

  return primary_particle;
}

std::set<PHG4Shower*> BaseTruthEval::all_secondary_showers(PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return std::set<PHG4Shower*>();
  }

  if (m_Strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++m_Errors;
    return std::set<PHG4Shower*>();
  }

  std::set<PHG4Shower*> subshowers;

  PHG4TruthInfoContainer::ShowerRange range = m_TruthInfo->GetSecondaryShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    PHG4Shower* shower = iter->second;

    if (m_Strict)
    {
      assert(shower);
    }
    else if (!shower)
    {
      ++m_Errors;
    }

    if (shower)
    {
      if (shower->get_parent_shower_id() == shower->get_id())
      {
        subshowers.insert(shower);
      }
    }
  }

  return subshowers;
}

bool BaseTruthEval::are_same_shower(PHG4Shower* s1, PHG4Shower* s2)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(s1);
    assert(s2);
  }
  else if (!s1 || !s2)
  {
    ++m_Errors;
    return false;
  }

  if (s1->get_id() == s2->get_id()) return true;
  return false;
}

bool BaseTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(p1);
    assert(p2);
  }
  else if (!p1 || !p2)
  {
    ++m_Errors;
    return false;
  }

  if (p1->get_track_id() == p2->get_track_id()) return true;
  return false;
}

bool BaseTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(vtx1);
    assert(vtx2);
  }
  else if (!vtx1 || !vtx2)
  {
    ++m_Errors;
    return false;
  }

  if (vtx1->get_id() == vtx2->get_id()) return true;
  return false;
}

PHG4Particle* BaseTruthEval::get_particle(PHG4Hit* g4hit)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++m_Errors;
    return nullptr;
  }

  PHG4Particle* particle = m_TruthInfo->GetParticle(g4hit->get_trkid());
  if (m_Strict)
  {
    assert(particle);
  }
  else if (!particle)
  {
    ++m_Errors;
  }

  return particle;
}

PHG4Shower* BaseTruthEval::get_primary_shower(PHG4Hit* g4hit)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++m_Errors;
    return nullptr;
  }

  PHG4Shower* shower = m_TruthInfo->GetShower(g4hit->get_shower_id());
  if (m_Strict)
  {
    assert(shower);
  }
  else if (!shower)
  {
    ++m_Errors;
  }

  return shower;
}

PHG4Particle* BaseTruthEval::get_primary_particle(PHG4Hit* g4hit)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return nullptr;
  }

  if (m_Strict)
  {
    assert(g4hit);
  }
  else if (!g4hit)
  {
    ++m_Errors;
    return nullptr;
  }

  PHG4Particle* particle = get_particle(g4hit);
  PHG4Particle* primary = get_primary_particle(particle);

  if (m_Strict)
  {
    assert(primary);
  }
  else if (!primary)
  {
    ++m_Errors;
  }

  return primary;
}

PHG4Particle* BaseTruthEval::get_particle(const int trackid)
{
  return m_TruthInfo->GetParticle(trackid);
}

bool BaseTruthEval::is_g4hit_from_primary_shower(PHG4Hit* g4hit, PHG4Shower* shower)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(g4hit);
    assert(shower);
  }
  else if (!g4hit || !shower)
  {
    ++m_Errors;
    return false;
  }

  if (g4hit->get_shower_id() == shower->get_id())
  {
    return true;
  }

  return false;
}

bool BaseTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle)
{
  if (!has_reduced_node_pointers())
  {
    ++m_Errors;
    return false;
  }

  if (m_Strict)
  {
    assert(g4hit);
    assert(particle);
  }
  else if (!g4hit || !particle)
  {
    ++m_Errors;
    return false;
  }

  if (g4hit->get_trkid() == particle->get_track_id())
  {
    return true;
  }

  return false;
}

void BaseTruthEval::get_node_pointers(PHCompositeNode* topNode)
{
  m_TruthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_TruthInfo)
  {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  return;
}

bool BaseTruthEval::has_reduced_node_pointers()
{
  if (m_Strict)
    assert(m_TruthInfo);
  else if (!m_TruthInfo)
    return false;

  return true;
}
