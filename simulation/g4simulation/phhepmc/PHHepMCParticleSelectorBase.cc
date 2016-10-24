#include "PHHepMCParticleSelectorBase.h"

//#include "PHG4Particlev1.h"

//#include "PHG4InEvent.h"
//#include "PHG4VtxPoint.h"
//#include "PHG4TruthInfoContainer.h"

#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

using namespace std;

PHHepMCParticleSelectorBase::PHHepMCParticleSelectorBase(const string &name):
  SubsysReco(name)
{
  return;
}

PHHepMCParticleSelectorBase::~PHHepMCParticleSelectorBase()
{
  return;
}

int
PHHepMCParticleSelectorBase::InitRun(PHCompositeNode *topNode)
{
  cout << PHWHERE << " " << Name() << " using empty InitRun()" << endl;
  return 0;
}

int PHHepMCParticleSelectorBase::process_event(PHCompositeNode *topNode)
{
  cout << PHWHERE << " " << Name() << " using empty process_event()" << endl;
  return 0;
}

void PHHepMCParticleSelectorBase::SetParticle(const int pid) 
{
  _theParticle = pid;
  return;
}

void PHHepMCParticleSelectorBase::AddParent(const int pid)
{
  _theParents.push_back(pid);
  return;
}

void PHHepMCParticleSelectorBase::AddDaughter(const int pid)
{
  _theDaughters.push_back(pid);
  return;
}

