#include "Fun4AllHepMCPileupInputManager.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"
#include "PHHepMCGenHelper.h"  // for PHHepMCGenHelper, PHHepMCGen...

#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY_SOME
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHRandomSeed.h>

#include <HepMC/GenEvent.h>
#include <HepMC/IO_GenEvent.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>  // for assert
#include <fstream>
#include <iostream>

Fun4AllHepMCPileupInputManager::Fun4AllHepMCPileupInputManager(
    const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllHepMCInputManager(name, nodename, topnodename)
{
  //! repeatedly read the input file
  Repeat(1);

  //! If set_embedding_id(i) with a negative number or 0, the pile up event will be inserted with increasing positive embedding_id. This is the default operation mode.
  PHHepMCGenHelper::set_embedding_id(-1);

  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(RandomGenerator, seed);

  return;
}

Fun4AllHepMCPileupInputManager::~Fun4AllHepMCPileupInputManager()
{
  gsl_rng_free(RandomGenerator);
}

int Fun4AllHepMCPileupInputManager::SkipForThisManager(const int nevents)
{
  for (int i = 0; i < nevents; ++i)
  {
    if (m_SignalInputManager)
    {
      m_SignalEventNumber = m_SignalInputManager->MyCurrentEvent(i);
    }
    int iret = run(1, true);
    if (iret)
    {
      return iret;
    }
  }
  return 0;
}

int Fun4AllHepMCPileupInputManager::run(const int /*nevents*/, const bool skip)
{
  if (_first_run)
  {
    _first_run = false;

    _ave_coll_per_crossing = _collision_rate * _time_between_crossings * 1e-9;
    _min_crossing = _min_integration_time / _time_between_crossings;
    _max_crossing = _max_integration_time / _time_between_crossings;

    if (Verbosity() >= VERBOSITY_SOME)
    {
      std::cout << "Fun4AllHepMCPileupInputManager::run:first event:  - ";
      std::cout << " _ave_coll_per_crossing = " << _ave_coll_per_crossing;
      std::cout << " _min_crossing = " << _min_crossing;
      std::cout << " _max_crossing = " << _max_crossing;
      std::cout << ". Start first event." << std::endl;
    }
  }
  if (m_SignalInputManager && !skip)
  {
    m_SignalEventNumber = m_SignalInputManager->MyCurrentEvent();
  }
  // if we have pushed back events we read them here and return
  if (m_EventPushedBackFlag)
  {
    // if m_EventPushedBackFlag = -1 we have no events but still created the temporary HepMC file
    // so we still need to remove it
    // if m_EventPushedBackFlag > 0 we have events pushed back which need to be recovered
    if (m_EventPushedBackFlag > 0)
    {
      HepMC::IO_GenEvent ascii_tmp_in(m_HepMCTmpFile, std::ios::in);
      HepMC::GenEvent *evttmp = ascii_tmp_in.read_next_event();
      while (ascii_tmp_in.rdstate() == 0)
      {
        if (evttmp->event_number() != m_SignalEventNumber)
        {
          double crossing_time = m_EventNumberMap.find(evttmp->event_number())->second;
          InsertEvent(evttmp, crossing_time);
        }
        evttmp = ascii_tmp_in.read_next_event();
      }
    }
    m_EventNumberMap.clear();
    m_EventPushedBackFlag = 0;
    remove(m_HepMCTmpFile.c_str());
    return 0;
  }
  // toss multiple crossings all the way back
  for (int icrossing = _min_crossing; icrossing <= _max_crossing; ++icrossing)
  {
    double crossing_time = _time_between_crossings * icrossing;

    int ncollisions = gsl_ran_poisson(RandomGenerator, _ave_coll_per_crossing);
    //    if (icrossing == 0) --ncollisions; // Jin: this is not necessary.
    // Triggering an interesting crossing does not change the distribution of the number of background collisions in that crossing

    for (int icollision = 0; icollision < ncollisions; ++icollision)
    {
      // loop until retrieve a valid event
      while (true)
      {
        if (!IsOpen())
        {
          if (FileListEmpty())
          {
            if (Verbosity() > 0)
            {
              std::cout << Name() << ": No Input file open" << std::endl;
            }
            return -1;
          }
          else
          {
            if (OpenNextFile())
            {
              std::cout << Name() << ": No Input file from filelist opened" << std::endl;
              return -1;
            }
          }
        }
        {
          if (ReadOscar())
          {
            evt = ConvertFromOscar();
            if (evt && m_SignalEventNumber == evt->event_number())
            {
              delete evt;
              evt = ConvertFromOscar();
            }
          }
          else
          {
            evt = ascii_in->read_next_event();
            if (evt && m_SignalEventNumber == evt->event_number())
            {
              delete evt;
              evt = ascii_in->read_next_event();
            }
          }
        }

        if (!evt)
        {
          if (Verbosity() > 1)
          {
            std::cout << "error type: " << ascii_in->error_type()
                 << ", rdstate: " << ascii_in->rdstate() << std::endl;
          }
          fileclose();
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << "Fun4AllHepMCPileupInputManager::run::" << Name();
            if (skip) std::cout << " skip";
            std::cout << " hepmc evt no: " << evt->event_number() << std::endl;
          }
          events_total++;
          events_thisfile++;
          // check if the local SubsysReco discards this event
          if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
          {
            delete evt;
            evt = nullptr;
            ResetEvent();
          }
          else
          {
            m_EventNumberMap.insert(std::make_pair(evt->event_number(), crossing_time));
            break;  // got the evt, move on
          }
        }
      }  // loop until retrieve a valid event
      if (!skip)
      {
        InsertEvent(evt, crossing_time);
      }
      else
      {
	delete evt;
	evt = nullptr;
      }
    }  //    for (int icollision = 0; icollision < ncollisions; ++icollision)

  }  //  for (int icrossing = _min_crossing; icrossing <= _max_crossing; ++icrossing)
  return 0;
}

int Fun4AllHepMCPileupInputManager::ResetEvent()
{
  m_EventNumberMap.clear();
  m_SignalEventNumber = 0;
  return 0;
}

int Fun4AllHepMCPileupInputManager::PushBackEvents(const int i)
{
  if (i == 1)
  {
    if (m_HepMCTmpFile.empty())
    {
      // we need to create this filename just once, we reuse it. Do it only if we need it
      m_HepMCTmpFile = "/tmp/HepMCTmpPileEvent-" + Name() + "-" + std::to_string(getpid()) + ".hepmc";
    }
    m_EventPushedBackFlag = -1;
    HepMC::IO_GenEvent ascii_io(m_HepMCTmpFile, std::ios::out);
    PHHepMCGenEventMap *geneventmap = PHHepMCGenHelper::get_geneventmap();
    for (auto iter = geneventmap->begin(); iter != geneventmap->end(); ++iter)
    {
      if (m_EventNumberMap.find((iter->second)->getEvent()->event_number()) != m_EventNumberMap.end())
      {
        m_EventPushedBackFlag = 1;
        HepMC::GenEvent *evttmp = (iter->second)->getEvent();
        ascii_io << evttmp;
      }
    }
    return 0;
  }
  return -1;
}
int Fun4AllHepMCPileupInputManager::InsertEvent(HepMC::GenEvent *evt, const double crossing_time)
{
  PHHepMCGenEventMap *geneventmap = PHHepMCGenHelper::get_geneventmap();
  PHHepMCGenEvent *genevent = nullptr;
  if (PHHepMCGenHelper::get_embedding_id() > 0)
  {
    //! If set_embedding_id(i) with a positive number, the pile up event will be inserted with increasing positive embedding_id. This would be a strange way to use pile up.

    genevent = geneventmap->insert_active_event(get_PHHepMCGenEvent_template() );
  }
  else
  {
    //! If set_embedding_id(i) with a negative number or 0, the pile up event will be inserted with increasing positive embedding_id. This is the default operation mode.

    genevent = geneventmap->insert_background_event(get_PHHepMCGenEvent_template() );
  }
  assert(genevent);
  assert(evt);
  genevent->addEvent(evt);
  PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation(genevent);
  // place to the crossing center in time
  genevent->moveVertex(0, 0, 0, crossing_time);
  return 0;
}
