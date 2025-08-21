//
// Inspired by code from ATLAS.  Thanks!
//
#include "HepMCFlowAfterBurner.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"

#include <flowafterburner/flowAfterburner.h> // Afterburner class
#include <flowafterburner/AfterburnerAlgo.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>               // for SubsysReco

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <CLHEP/Random/MTwistEngine.h>
#include <CLHEP/Random/RandomEngine.h>

#include <iostream>
#include <iterator>                           // for operator!=, reverse_ite...
#include <set>                                // for set, _Rb_tree_const_ite...
#include <string>
#include <utility>                            // for pair


using namespace std;

class PHCompositeNode;
namespace HepMC { class GenEvent; }



// we want to keep the default eta range identical between here and 
// the flowAfterburner executable. If you change the default eta range here
// please apply the same change to generators/flowAfterburner/main.cc
HepMCFlowAfterBurner::HepMCFlowAfterBurner(const std::string &name)
  : SubsysReco(name)
{
}

int HepMCFlowAfterBurner::Init(PHCompositeNode */*topNode*/)
{
  if (seedset)
  {
    randomSeed = seed;
  }
  else
  {
    randomSeed = PHRandomSeed();
  }

  m_engine = new CLHEP::MTwistEngine(randomSeed);
  m_afterburner = new Afterburner(algorithmName, m_engine, mineta, maxeta, minpt, maxpt);
  m_flowalgo = m_afterburner->getAlgo();
  // you can set other algo parameters here if needed
  if (enableFlucuations)
  {
    m_flowalgo->enable_fluctuations();
  }
  
  for (unsigned int i = 1; i <= 6; ++i)
  {//apply the scale factors to the flow harmonics if any
    if (flowScales[i - 1] != 1.0F)
    {
      m_flowalgo->set_single_scale_N(i, flowScales[i - 1]);
    }  
  }

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


    m_afterburner->flowAfterburner(evt);

    for ( unsigned int i=1; i<=6; ++i)
    {
      genevt->set_flow_psi(i, m_afterburner->getPsiN(i));
      if (Verbosity() > 1)
      {
        cout << "  set reaction plane angle psi_" << i << " = " << genevt->get_flow_psi(i) << endl;
      }
    }

    if (Verbosity() > 1)
    {
      m_afterburner->getAlgo()->print();
    }

    // flowAfterburner(evt, engine, algorithmName, mineta, maxeta, minpt, maxpt);
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
  if (m_engine)
  {
    delete m_engine;
    m_engine = new CLHEP::MTwistEngine(randomSeed);
  }
  return;
}

void HepMCFlowAfterBurner::SaveRandomState(const string &savefile)
{
  if (m_engine)
  {
    m_engine->saveStatus(savefile.c_str());
    return;
  }
  cout << PHWHERE << " Random engine not started yet" << endl;
}

void HepMCFlowAfterBurner::RestoreRandomState(const string &savefile)
{
  if (m_engine)
  {
    m_engine->restoreStatus(savefile.c_str());
    return;
  }
  cout << PHWHERE << " Random engine not started yet" << endl;
}

void HepMCFlowAfterBurner::Print(const string &/*what*/) const
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

  m_flowalgorithm = AfterburnerAlgo::getAlgoFromName(name); // will print error message if algo name is unknown
  algorithmName = AfterburnerAlgo::getAlgoName(m_flowalgorithm); // make sure the name is consistent
  return;
}

void HepMCFlowAfterBurner::scaleFlow(const float scale, const unsigned int n)
{
  if ( n == 0 )
  { // set all scales
    for (unsigned int i = 0; i < 6; ++i)
    {
      flowScales[i] = scale;
    }
  }
  else if ( n > 0 && n <= 6 )
  { // set specific harmonic
    flowScales[n - 1] = scale;
  }
  else
  { // out of range
    std::cout << "HepMCFlowAfterBurner::scaleFlow - ERROR: n = " << n << " is out of range.  Must be between 0 (all) or 1,..,6." << std::endl;
  }
}

