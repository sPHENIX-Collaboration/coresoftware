#include "FermimotionAfterburner.h"

#include "FermiMotion.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>               // for SubsysReco

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <gsl/gsl_rng.h>

#include <iostream>
#include <set>                                // for set, _Rb_tree_const_ite...
#include <string>
#include <utility>   

namespace HepMC { class GenEvent; }

//____________________________________________________________________________..
FermimotionAfterburner::FermimotionAfterburner(const std::string &name):
  SubsysReco(name)
  
{
  
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

//____________________________________________________________________________..
FermimotionAfterburner::~FermimotionAfterburner()
{
  gsl_rng_free(RandomGenerator);
}

//____________________________________________________________________________..
int FermimotionAfterburner::Init(PHCompositeNode *topNode)
{
  
  unsigned int seed = PHRandomSeed();  
  gsl_rng_set(RandomGenerator, seed);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FermimotionAfterburner::InitRun(PHCompositeNode *topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FermimotionAfterburner::process_event(PHCompositeNode *topNode)
{
  std::cout << "process_event(PHCompositeNode *topNode) Processing Event Shuhang testing" << std::endl;
  AddpF(topNode);
  
   
  return Fun4AllReturnCodes::EVENT_OK;
  
  
}

//____________________________________________________________________________..
int FermimotionAfterburner::ResetEvent(PHCompositeNode *topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FermimotionAfterburner::EndRun(const int runnumber)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FermimotionAfterburner::End(PHCompositeNode *topNode)
{
  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FermimotionAfterburner::Reset(PHCompositeNode *topNode)
{
  
 
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void FermimotionAfterburner::Print(const std::string &what) const
{
 
}

//____________________________________________________________________________..


void FermimotionAfterburner::AddpF(PHCompositeNode *topNode)

{
PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
 std::cout<<"looping over events"<<std::endl;
 for (PHHepMCGenEventMap::Iter iter = genevtmap->begin(); iter != genevtmap->end(); ++iter)
    {
      std::cout<<"getting event"<<std::endl;
      PHHepMCGenEvent *genevt = iter->second;
      HepMC::GenEvent *evt = genevt->getEvent();
      std::cout<<"done"<<std::endl;
      if (!evt)
	{
	  std::cout << PHWHERE << " no evt pointer under HEPMC Node found" << std::endl;
	}
      std::cout<<"applying fermimotion"<<std::endl;
      FermiMotion(evt,RandomGenerator);
   
    }
  
 
}
