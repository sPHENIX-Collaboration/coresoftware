#include "FermimotionAfterburner.h"
#include "FermiMotion.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>               // for SubsysReco

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <iostream>
#include <iterator>                           // for operator!=, reverse_ite...
#include <set>                                // for set, _Rb_tree_const_ite...
#include <string>
#include <utility>   

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h> 

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>  // for GenParticle
#include <HepMC/GenRanges.h>
#include <HepMC/GenVertex.h>      // for GenVertex, GenVertex::part...
#include <HepMC/HeavyIon.h>       // for HeavyIon
#include <HepMC/IteratorRange.h>  // for children, descendants
#include <HepMC/SimpleVector.h>   // for FourVector

#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Vector/LorentzVector.h>

#include <cmath>
#include <map>  // for map

using namespace std;
using namespace HepMC;

namespace HepMC { class GenEvent; }
CLHEP::HepRandomEngine *engine = nullptr;
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
  cout << "process_event(PHCompositeNode *topNode) Processing Event Shuhang testing" << endl;
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
 cout<<"looping over events"<<endl; 
 for (PHHepMCGenEventMap::Iter iter = genevtmap->begin(); iter != genevtmap->end(); ++iter)
    {
      cout<<"getting event"<<endl;
      PHHepMCGenEvent *genevt = iter->second;
      HepMC::GenEvent *evt = genevt->getEvent();
      cout<<"done"<<endl;
      if (!evt)
	{
	  cout << PHWHERE << " no evt pointer under HEPMC Node found" << endl;
	}
      cout<<"applying fermimotion"<<std::endl;
      FermiMotion(evt,RandomGenerator);
   
    }
  
 
}
