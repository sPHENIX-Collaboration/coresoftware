#include "PHPythia8.h"

#include "PHPy8GenTrigger.h"

#include <phhepmc/PHGenIntegral.h>  // for PHGenIntegral
#include <phhepmc/PHGenIntegralv1.h>
#include <phhepmc/PHHepMCGenHelper.h>  // for PHHepMCGenHelper

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

#include <HepMC/GenEvent.h>
#include <HepMC/Units.h>            // for GEV, MM
#include <HepMC/WeightContainer.h>  // for WeightContainer

#include <Pythia8/Event.h>  // for Event
#include <Pythia8/Info.h>   // for Info
#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>

#include <cassert>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>  // for operator<<, endl

/**
 * @brief Construct a PHPythia8 generator instance and configure HepMC conversion.
 *
 * Initializes the Pythia8 engine using the path from the environment variable
 * `PYTHIA8`, configures a HepMC::Pythia8ToHepMC converter to store process,
 * PDF, and cross-section information, and sets the default embedding ID to 1.
 * The constructor preserves and restores std::cout formatting around Pythia8
 * construction to avoid altering global stream state.
 *
 * If `PYTHIA8` is not set, an error message is printed and the Pythia8 instance
 * remains uninitialized.
 *
 * @param name Name forwarded to the SubsysReco base class (module instance name).
 */
PHPythia8::PHPythia8(const std::string &name)
  : SubsysReco(name)
{
  char *charPath = getenv("PYTHIA8");
  if (!charPath)
  {
    std::cout << "PHPythia8::Could not find $PYTHIA8 path!" << std::endl;
    return;
  }

  std::string thePath(charPath);
  thePath += "/xmldoc/";
  // the pythia8 ctor messes with the formatting, so we save the cout state here
  // and restore it later
  std::ios old_state(nullptr);
  old_state.copyfmt(std::cout);
  m_Pythia8.reset(new Pythia8::Pythia(thePath));
  std::cout.copyfmt(old_state);
  m_Pythia8ToHepMC.reset(new HepMC::Pythia8ToHepMC());
  m_Pythia8ToHepMC->set_store_proc(true);
  m_Pythia8ToHepMC->set_store_pdf(true);
  m_Pythia8ToHepMC->set_store_xsec(true);

  PHHepMCGenHelper::set_embedding_id(1);  // default embedding ID to 1
}

/**
 * @brief Initialize the Pythia8 generator, configure nodes, and seed the RNG.
 *
 * Performs module initialization: reads an optional configuration file and any
 * queued Pythia command strings, creates the required node tree under the
 * provided top-level node, sets Pythia's random seed (mapped from PHRandomSeed
 * into Pythia's valid range) and prints it for reproducibility, then calls
 * Pythia8::init().
 *
 * @param topNode Top-level PHCompositeNode under which generator nodes are created.
 * @return int Fun4All return code; returns Fun4AllReturnCodes::EVENT_OK on success.
 */
int PHPythia8::Init(PHCompositeNode *topNode)
{
  if (!m_ConfigFileName.empty())
  {
    read_config(m_ConfigFileName);
  }
  for (auto &m_Command : m_Commands)
  {
    m_Pythia8->readString(m_Command);
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
    m_Pythia8->readString(std::format("Random:seed = {}", seed));
  }
  else
  {
    std::cout << PHWHERE << " ERROR: seed " << seed << " is not valid" << std::endl;
    exit(1);
  }
  // print out seed so we can make this is reproducible
  std::cout << "PHPythia8 random seed: " << seed << std::endl;


  // this is empirical - something in the pythia8::init() method interferes
  // with our macros (it sets the tpc drift verlocity back to 0)
  // not the feintest idea right now what this could be
  // but saving the old cout state and restoring it aftwerwards
  // gets our tpc drift velocity back
  std::ios old_state(nullptr);
  old_state.copyfmt(std::cout);

  m_Pythia8->init();

  std::cout.copyfmt(old_state);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    std::cout << "PHPythia8::End - I'm here!" << std::endl;
  }

  if (Verbosity() >= VERBOSITY_SOME)
  {
    //-* dump out closing info (cross-sections, etc)
    m_Pythia8->stat();

    // match pythia printout
    std::cout << " |                                                                "
              << "                                                 | " << std::endl;
    std::cout << "                         PHPythia8::End - " << m_EventCount
              << " events passed trigger" << std::endl;
    std::cout << "                         Fraction passed: " << m_EventCount
              << "/" << m_Pythia8->info.nAccepted()
              << " = " << m_EventCount / float(m_Pythia8->info.nAccepted()) << std::endl;
    std::cout << " *-------  End PYTHIA Trigger Statistics  ------------------------"
              << "-------------------------------------------------* " << std::endl;

    if (m_IntegralNode)
    {
      std::cout << "Integral information on stored on node RUN/PHGenIntegral:" << std::endl;
      m_IntegralNode->identify();
      std::cout << " *-------  End PYTHIA Integral Node Print  ------------------------"
                << "-------------------------------------------------* " << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________________________
int PHPythia8::read_config(const std::string &cfg_file)
{
  m_ConfigFileName = cfg_file;

  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "PHPythia8::read_config - Reading " << m_ConfigFileName << std::endl;
  }

  std::ifstream infile(m_ConfigFileName);
  if (infile.fail())
  {
    std::cout << "PHPythia8::read_config - Failed to open file " << m_ConfigFileName << std::endl;
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
    std::cout << "PHPythia8::process_event - event: " << m_EventCount << std::endl;
  }

  bool passedGen = false;
  bool passedTrigger = false;
  //  int genCounter = 0;

  while (!passedTrigger)
  {
    //    ++genCounter;

    // generate another pythia event
    while (!passedGen)
    {
      passedGen = m_Pythia8->next();
    }

    // test trigger logic
    bool andScoreKeeper = true;
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      std::cout << "PHPythia8::process_event - triggersize: " << m_RegisteredTriggers.size() << std::endl;
    }

    for (auto &m_RegisteredTrigger : m_RegisteredTriggers)
    {
      bool trigResult = m_RegisteredTrigger->Apply(m_Pythia8.get());

      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "PHPythia8::process_event trigger: "
                  << m_RegisteredTrigger->GetName() << "  " << trigResult << std::endl;
      }

      if (m_TriggersOR && trigResult)
      {
        passedTrigger = true;
        break;
      }
      if (m_TriggersAND)
      {
        andScoreKeeper &= trigResult;
      }

      if (Verbosity() >= VERBOSITY_EVEN_MORE && !passedTrigger)
      {
        std::cout << "PHPythia8::process_event - failed trigger: "
                  << m_RegisteredTrigger->GetName() << std::endl;
      }
    }

    if ((andScoreKeeper && m_TriggersAND) || (m_RegisteredTriggers.empty()))
    {
      passedTrigger = true;
      //      genCounter = 0;
    }

    passedGen = false;
  }

  // print
  if (Verbosity())
  {
    m_Pythia8->event.list();
  }

  // fill HepMC object with event & pass to

  auto *genevent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  m_Pythia8ToHepMC->fill_next_event(*m_Pythia8, genevent, m_EventCount);
  // Enable continuous reweighting by storing additional reweighting factor
  if (m_SaveEventWeightFlag)
  {
    genevent->weights().push_back(m_Pythia8->info.weight());
  }

  /* pass HepMC to PHNode*/
  auto *success = PHHepMCGenHelper::insert_event(genevent);
  if (!success)
  {
    std::cout << "PHPythia8::process_event - Failed to add event to HepMC record!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // print outs

  if (Verbosity() >= VERBOSITY_MORE)
  {
    std::cout << "PHPythia8::process_event - FINISHED WHOLE EVENT" << std::endl;
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
    std::cout << PHWHERE << "RUN Node missing doing nothing" << std::endl;
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
      std::cout << "PHPythia8::create_node_tree - Fatal Error - "
                << "RUN/PHGenIntegral node already exist. "
                << "It is messy to overwrite integrated luminosities. Please turn off this function in the macro with " << std::endl;
      std::cout << "                              PHPythia8::save_integrated_luminosity(false);" << std::endl;
      std::cout << "The current RUN/PHGenIntegral node is ";
      m_IntegralNode->identify(std::cout);

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
    std::cout << "PHPythia8::registerTrigger - trigger " << theTrigger->GetName() << " registered" << std::endl;
  }
  m_RegisteredTriggers.push_back(theTrigger);
}