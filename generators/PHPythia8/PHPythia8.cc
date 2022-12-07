#include "PHPythia8.h"

#include "PHPy8GenTrigger.h"

#include <phhepmc/PHGenIntegral.h>  // for PHGenIntegral
#include <phhepmc/PHGenIntegralv1.h>
#include <phhepmc/PHHepMCGenHelper.h>  // for PHHepMCGenHelper

#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBO...
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop

#include <HepMC/Units.h>            // for GEV, MM
#include <HepMC/WeightContainer.h>  // for WeightContainer

#include <Pythia8/Event.h>  // for Event
#include <Pythia8/Info.h>   // for Info
#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>

#include <boost/format.hpp>

#include <cassert>
#include <cstdlib>
#include <iostream>  // for operator<<, endl

class PHHepMCGenEvent;

using namespace std;

PHPythia8::PHPythia8(const std::string &name)
  : SubsysReco(name)
  , m_EventCount(0)
  , m_TriggersOR(true)
  , m_TriggersAND(false)
  , m_Pythia8(nullptr)
  , m_ConfigFileName("phpythia8.cfg")
  , m_Pythia8ToHepMC(nullptr)
  , m_SaveEventWeightFlag(true)
  , m_SaveIntegratedLuminosityFlag(true)
  , m_IntegralNode(nullptr)
{
  char *charPath = getenv("PYTHIA8");
  if (!charPath)
  {
    cout << "PHPythia8::Could not find $PYTHIA8 path!" << endl;
    return;
  }

  std::string thePath(charPath);
  thePath += "/xmldoc/";
  m_Pythia8 = new Pythia8::Pythia(thePath.c_str());

  m_Pythia8ToHepMC = new HepMC::Pythia8ToHepMC();
  m_Pythia8ToHepMC->set_store_proc(true);
  m_Pythia8ToHepMC->set_store_pdf(true);
  m_Pythia8ToHepMC->set_store_xsec(true);

  PHHepMCGenHelper::set_embedding_id(1);  // default embedding ID to 1
}

PHPythia8::~PHPythia8()
{
  delete m_Pythia8;
  delete m_Pythia8ToHepMC;
}

int PHPythia8::Init(PHCompositeNode *topNode)
{
  if (!m_ConfigFileName.empty())
  {
    read_config(m_ConfigFileName);
  }
  for (unsigned int j = 0; j < m_Commands.size(); j++)
  {
    m_Pythia8->readString(m_Commands[j]);
  }

  create_node_tree(topNode);

  // PYTHIA8 has very specific requires for its random number range
  // I map the designated unique seed from recoconst into something
  // acceptable for PYTHIA8

  unsigned int seed = PHRandomSeed();

  if (seed > 900000000)
  {
    seed = seed % 900000000;
  }

  if ((seed > 0) && (seed <= 900000000))
  {
    m_Pythia8->readString("Random:setSeed = on");
    m_Pythia8->readString(str(boost::format("Random:seed = %1%") % seed));
  }
  else
  {
    cout << PHWHERE << " ERROR: seed " << seed << " is not valid" << endl;
    exit(1);
  }
  // print out seed so we can make this is reproducible
  cout << "PHPythia8 random seed: " << seed << endl;

  m_Pythia8->init();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    cout << "PHPythia8::End - I'm here!" << endl;
  }

  if (Verbosity() >= VERBOSITY_SOME)
  {
    //-* dump out closing info (cross-sections, etc)
    m_Pythia8->stat();

    //match pythia printout
    cout << " |                                                                "
         << "                                                 | " << endl;
    cout << "                         PHPythia8::End - " << m_EventCount
         << " events passed trigger" << endl;
    cout << "                         Fraction passed: " << m_EventCount
         << "/" << m_Pythia8->info.nAccepted()
         << " = " << m_EventCount / float(m_Pythia8->info.nAccepted()) << endl;
    cout << " *-------  End PYTHIA Trigger Statistics  ------------------------"
         << "-------------------------------------------------* " << endl;

    if (m_IntegralNode)
    {
      cout << "Integral information on stored on node RUN/PHGenIntegral:" << endl;
      m_IntegralNode->identify();
      cout << " *-------  End PYTHIA Integral Node Print  ------------------------"
           << "-------------------------------------------------* " << endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________________________
int PHPythia8::read_config(const string &cfg_file)
{
  m_ConfigFileName = cfg_file;

  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "PHPythia8::read_config - Reading " << m_ConfigFileName << endl;
  }

  ifstream infile(m_ConfigFileName);
  if (infile.fail())
  {
    cout << "PHPythia8::read_config - Failed to open file " << m_ConfigFileName << endl;
    exit(2);
  }

  m_Pythia8->readFile(m_ConfigFileName);

  return Fun4AllReturnCodes::EVENT_OK;
}

//-* print pythia config info
void PHPythia8::print_config() const
{
  m_Pythia8->info.list();
}

int PHPythia8::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    cout << "PHPythia8::process_event - event: " << m_EventCount << endl;
  }

  bool passedGen = false;
  bool passedTrigger = false;
  int genCounter = 0;

  while (!passedTrigger)
  {
    ++genCounter;

    // generate another pythia event
    while (!passedGen)
    {
      passedGen = m_Pythia8->next();
    }

    // test trigger logic

    bool andScoreKeeper = true;
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      cout << "PHPythia8::process_event - triggersize: " << m_RegisteredTriggers.size() << endl;
    }

    for (unsigned int tr = 0; tr < m_RegisteredTriggers.size(); tr++)
    {
      bool trigResult = m_RegisteredTriggers[tr]->Apply(m_Pythia8);

      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "PHPythia8::process_event trigger: "
             << m_RegisteredTriggers[tr]->GetName() << "  " << trigResult << endl;
      }

      if (m_TriggersOR && trigResult)
      {
        passedTrigger = true;
        break;
      }
      else if (m_TriggersAND)
      {
        andScoreKeeper &= trigResult;
      }

      if (Verbosity() >= VERBOSITY_EVEN_MORE && !passedTrigger)
      {
        cout << "PHPythia8::process_event - failed trigger: "
             << m_RegisteredTriggers[tr]->GetName() << endl;
      }
    }

    if ((andScoreKeeper && m_TriggersAND) || (m_RegisteredTriggers.size() == 0))
    {
      passedTrigger = true;
      genCounter = 0;
    }

    passedGen = false;
  }

  // fill HepMC object with event & pass to

  HepMC::GenEvent *genevent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  m_Pythia8ToHepMC->fill_next_event(*m_Pythia8, genevent, m_EventCount);
  // Enable continuous reweighting by storing additional reweighting factor
  if (m_SaveEventWeightFlag)
  {
    genevent->weights().push_back(m_Pythia8->info.weight());
  }

  /* pass HepMC to PHNode*/
  PHHepMCGenEvent *success = PHHepMCGenHelper::insert_event(genevent);
  if (!success)
  {
    cout << "PHPythia8::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // print outs

  if (Verbosity() >= VERBOSITY_MORE)
  {
    cout << "PHPythia8::process_event - FINISHED WHOLE EVENT" << endl;
  }
  if (m_EventCount < 2 && Verbosity() >= VERBOSITY_SOME)
  {
    m_Pythia8->event.list();
  }
  if (m_EventCount >= 2 && Verbosity() >= VERBOSITY_A_LOT)
  {
    m_Pythia8->event.list();
  }

  ++m_EventCount;

  // save statistics
  if (m_IntegralNode)
  {
    m_IntegralNode->set_N_Generator_Accepted_Event(m_Pythia8->info.nAccepted());
    m_IntegralNode->set_N_Processed_Event(m_EventCount);
    m_IntegralNode->set_Sum_Of_Weight(m_Pythia8->info.weightSum());
    m_IntegralNode->set_Integrated_Lumi(m_Pythia8->info.nAccepted() / (m_Pythia8->info.sigmaGen() * 1e9));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::create_node_tree(PHCompositeNode *topNode)
{
  // HepMC IO
  PHHepMCGenHelper::create_node_tree(topNode);

  PHNodeIterator iter(topNode);
  PHCompositeNode *sumNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!sumNode)
  {
    cout << PHWHERE << "RUN Node missing doing nothing" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (m_SaveIntegratedLuminosityFlag)
  {
    m_IntegralNode = findNode::getClass<PHGenIntegral>(sumNode, "PHGenIntegral");
    if (!m_IntegralNode)
    {
      m_IntegralNode = new PHGenIntegralv1("PHPythia8 with embedding ID of " + std::to_string(PHHepMCGenHelper::get_embedding_id()));
      PHIODataNode<PHObject> *newmapnode = new PHIODataNode<PHObject>(m_IntegralNode, "PHGenIntegral", "PHObject");
      sumNode->addNode(newmapnode);
    }
    else
    {
      cout << "PHPythia8::create_node_tree - Fatal Error - "
           << "RUN/PHGenIntegral node already exist. "
           << "It is messy to overwrite integrated luminosities. Please turn off this function in the macro with " << endl;
      cout << "                              PHPythia8::save_integrated_luminosity(false);" << endl;
      cout << "The current RUN/PHGenIntegral node is ";
      m_IntegralNode->identify(cout);

      exit(EXIT_FAILURE);
    }
    assert(m_IntegralNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHPythia8::register_trigger(PHPy8GenTrigger *theTrigger)
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    cout << "PHPythia8::registerTrigger - trigger " << theTrigger->GetName() << " registered" << endl;
  }
  m_RegisteredTriggers.push_back(theTrigger);
}
