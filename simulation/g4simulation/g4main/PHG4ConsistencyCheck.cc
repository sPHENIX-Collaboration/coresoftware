#include "PHG4ConsistencyCheck.h"

#include "PHG4Hit.h"
#include "PHG4HitContainer.h"
#include "PHG4Particle.h"
#include "PHG4TruthInfoContainer.h"

#include <phool/getClass.h>

#include <iostream>  // for operator<<, basic_ostream::opera...
#include <map>       // for _Rb_tree_const_iterator, map<>::...
#include <set>       // for set
#include <utility>   // for pair

class PHCompositeNode;

PHG4ConsistencyCheck::PHG4ConsistencyCheck(const std::string &name)
  : SubsysReco(name)
  , errorcnt(0)
{
}

int PHG4ConsistencyCheck::process_event(PHCompositeNode *topNode)
{
  PHG4TruthInfoContainer *truthcont = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthcont)
  {
    return 0;
  }
  PHG4TruthInfoContainer::ConstRange trange = truthcont->GetParticleRange();
  PHG4TruthInfoContainer::ConstIterator titer;
  int imax = 1000000;
  for (titer = trange.first; titer != trange.second; ++titer)
  {
    if (titer->first < imax)
    {
      imax = titer->first;
    }
  }
  std::cout << "min index: " << imax << std::endl;
  std::pair<std::map<int, int>::const_iterator, std::map<int, int>::const_iterator> embtrk_b_e = truthcont->GetEmbeddedTrkIds();
  std::map<int, int>::const_iterator embiter;
  for (embiter = embtrk_b_e.first; embiter != embtrk_b_e.second; ++embiter)
  {
    std::cout << "embedded trkid: " << embiter->first << std::endl;
  }
  PHG4HitContainer *ghit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_CEMC_E");
  if (ghit)
  {
    PHG4HitContainer::ConstIterator hit;
    PHG4HitContainer::ConstRange hit_begin_end = ghit->getHits();
    std::set<int> printpart;
    for (hit = hit_begin_end.first; hit != hit_begin_end.second; ++hit)
    {
      int trkid = hit->second->get_trkid();
      PHG4Particle *part = truthcont->GetParticle(trkid);
      if (!part)
      {
        hit->second->identify();
        std::cout << "could not locate geant particle " << trkid << " in G4HIT_CEMC_E" << std::endl;
        errorcnt++;
      }
      else
      {
        int primary_id = part->get_primary_id();
        if (truthcont->isEmbeded(primary_id) > 0)
        {
          if (printpart.find(primary_id) == printpart.end())
          {
            std::cout << "primary id " << primary_id << " is embedded" << std::endl;
            printpart.insert(primary_id);
            PHG4Particle *parta = truthcont->GetParticle(primary_id);
            parta->identify();
          }
        }
      }
    }
  }
  ghit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_SVTX");
  if (ghit)
  {
    PHG4HitContainer::ConstIterator hit;
    PHG4HitContainer::ConstRange hit_begin_end = ghit->getHits();
    for (hit = hit_begin_end.first; hit != hit_begin_end.second; ++hit)
    {
      int trkid = hit->second->get_trkid();
      PHG4Particle *part = truthcont->GetParticle(trkid);
      if (!part)
      {
        std::cout << "could not locate geant particle " << trkid << " in G4HIT_SVTX" << std::endl;
        errorcnt++;
      }
    }
  }

  return 0;
}
