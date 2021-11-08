#include "FermimotionAfterburner.h"

#include "FermiMotion.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <gsl/gsl_rng.h>

#include <iostream>
#include <set>  // for set, _Rb_tree_const_ite...
#include <string>
#include <utility>

namespace HepMC
{
  class GenEvent;
}

//____________________________________________________________________________..
FermimotionAfterburner::FermimotionAfterburner(const std::string &name)
  : SubsysReco(name)

{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

//____________________________________________________________________________..
FermimotionAfterburner::~FermimotionAfterburner()
{
  gsl_rng_free(RandomGenerator);
}

//____________________________________________________________________________..
int FermimotionAfterburner::Init(PHCompositeNode */*topNode*/)
{
  unsigned int seed = PHRandomSeed();
  gsl_rng_set(RandomGenerator, seed);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FermimotionAfterburner::process_event(PHCompositeNode *topNode)
{
  AddpF(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void FermimotionAfterburner::AddpF(PHCompositeNode *topNode)
{
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  for (PHHepMCGenEventMap::Iter iter = genevtmap->begin(); iter != genevtmap->end(); ++iter)
  {
    PHHepMCGenEvent *genevt = iter->second;
    HepMC::GenEvent *evt = genevt->getEvent();
    if (!evt)
    {
      std::cout << PHWHERE << " no evt pointer under HEPMC Node found" << std::endl;
    }
    FermiMotion(evt, RandomGenerator);
  }
}
