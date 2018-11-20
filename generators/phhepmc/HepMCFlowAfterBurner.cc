//
// Inspired by code from ATLAS.  Thanks!
//
#include "HepMCFlowAfterBurner.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"

#include <flowafterburner/flowAfterburner.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <HepMC/GenEvent.h>

#include <CLHEP/Random/MTwistEngine.h>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandomEngine.h>

#include <iostream>
#include <string>

using namespace std;

CLHEP::HepRandomEngine *engine = nullptr;

set<string> algoset = {"MINBIAS", "MINBIAS_V2_ONLY", "CUSTOM"};

HepMCFlowAfterBurner::HepMCFlowAfterBurner(const std::string &name)
  : SubsysReco(name)
  , algorithmName("MINBIAS")
  , mineta(-2)
  , maxeta(2)
  , minpt(0.)
  , maxpt(100.)
  , seedset(0)
  , seed(0)
  , randomSeed(11793)
{
}

int HepMCFlowAfterBurner::Init(PHCompositeNode *topNode)
{
  if (seedset)
  {
    randomSeed = seed;
  }
  else
  {
    randomSeed = PHRandomSeed();
  }

  engine = new CLHEP::MTwistEngine(randomSeed);

  return 0;
}

int HepMCFlowAfterBurner::process_event(PHCompositeNode *topNode)
{
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  for (PHHepMCGenEventMap::ReverseIter iter = genevtmap->rbegin(); iter != genevtmap->rend(); ++iter)
  {
    PHHepMCGenEvent *genevt = iter->second;

    HepMC::GenEvent *evt = genevt->getEvent();
    if (!evt)
    {
      cout << PHWHERE << " no evt pointer under HEPMC Node found" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if (Verbosity() > 0)
    {
      cout << "calling flowAfterburner with algorithm "
           << algorithmName << ", mineta " << mineta
           << ", maxeta: " << maxeta << ", minpt: " << minpt
           << ", maxpt: " << maxpt << endl;
    }
    flowAfterburner(evt, engine, algorithmName, mineta, maxeta, minpt, maxpt);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void HepMCFlowAfterBurner::setSeed(const long il)
{
  seedset = 1;
  seed = il;
  randomSeed = seed;
  // just in case we are already running, kill the engine and make
  // a new one using the selected seed
  if (engine)
  {
    delete engine;
    engine = new CLHEP::MTwistEngine(randomSeed);
  }
  return;
}

void HepMCFlowAfterBurner::SaveRandomState(const string &savefile)
{
  if (engine)
  {
    engine->saveStatus(savefile.c_str());
    return;
  }
  cout << PHWHERE << " Random engine not started yet" << endl;
}

void HepMCFlowAfterBurner::RestoreRandomState(const string &savefile)
{
  if (engine)
  {
    engine->restoreStatus(savefile.c_str());
    return;
  }
  cout << PHWHERE << " Random engine not started yet" << endl;
}

void HepMCFlowAfterBurner::Print(const string &what) const
{
  cout << "FlowAfterBurner parameters:" << endl;
  cout << "algorithm: " << algorithmName << endl;
  cout << "mineta: " << mineta << ", maxeta: " << maxeta << endl;
  cout << "minpt: " << minpt << ", maxpt: " << maxpt << endl;
  cout << "Implemented algorithms: MINBIAS (default), MINBIAS_V2_ONLY, CUSTOM"
       << endl;
  return;
}

void HepMCFlowAfterBurner::setAlgorithmName(const std::string &name)
{
  auto it = algoset.find(name);
  if (it != algoset.end())
  {
    algorithmName = *it;
  }
  else
  {
    cout << "algorithm " << name << " not in list of possible algorithms" << endl;
      cout << "possible algorithms are" << endl;
      for (auto &al : algoset)
      {
	cout << al << endl;
      }
  }
  return;
}
