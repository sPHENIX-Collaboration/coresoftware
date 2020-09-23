#include "Fun4AllHepMCPileupInputManager.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"
#include "PHHepMCGenHelper.h"            // for PHHepMCGenHelper, PHHepMCGen...

#include <fun4all/Fun4AllBase.h>         // for Fun4AllBase::VERBOSITY_SOME
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <phool/PHRandomSeed.h>

#include <HepMC/GenEvent.h>
#include <HepMC/IO_GenEvent.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>                      // for assert
#include <fstream>
#include <iostream>

using namespace std;

Fun4AllHepMCPileupInputManager::Fun4AllHepMCPileupInputManager(
    const string &name, const string &nodename, const string &topnodename)
  : Fun4AllHepMCInputManager(name, nodename, topnodename)
  , _min_integration_time(-17500.0)
  , _max_integration_time(+17500.0)
  , _collision_rate(100.0e3)
  , _time_between_crossings(106.0)
  , _ave_coll_per_crossing(1.0)
  ,  // recalculated
  _min_crossing(0)
  ,  // recalculated
  _max_crossing(0)
  ,  // recalculated
  _first_run(true)
{
  //! repeatedly read the input file
  Repeat(1);

  //! If set_embedding_id(i) with a negative number or 0, the pile up event will be inserted with increasing positive embedding_id. This is the default operation mode.
  hepmc_helper.set_embedding_id(-1);

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
  for (int i=0; i<nevents; ++i)
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

int Fun4AllHepMCPileupInputManager::run(const int nevents, const bool skip)
{
  if (_first_run)
  {
    _first_run = false;

    _ave_coll_per_crossing = _collision_rate * _time_between_crossings * 1e-9;
    _min_crossing = _min_integration_time / _time_between_crossings;
    _max_crossing = _max_integration_time / _time_between_crossings;

    if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "Fun4AllHepMCPileupInputManager::run:first event:  - ";
      cout << " _ave_coll_per_crossing = " << _ave_coll_per_crossing;
      cout << " _min_crossing = " << _min_crossing;
      cout << " _max_crossing = " << _max_crossing;
      cout << ". Start first event." << endl;
    }
  }
  if (m_SignalInputManager && !skip)
  {
    m_SignalEventNumber = m_SignalInputManager->MyCurrentEvent();
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
      double t0 = crossing_time;

      // loop until retrieve a valid event
      while (true)
      {
        if (!IsOpen())
        {
          if (FileListEmpty())
          {
            if (Verbosity() > 0)
            {
              cout << Name() << ": No Input file open" << endl;
            }
            return -1;
          }
          else
          {
            if (OpenNextFile())
            {
              cout << Name() << ": No Input file from filelist opened" << endl;
              return -1;
            }
          }
        }

        if (save_evt)
        {  // if an event was pushed back, copy saved pointer and
           // reset save_evt pointer
          evt = save_evt;
          save_evt = nullptr;
        }
        else
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
            cout << "error type: " << ascii_in->error_type()
                 << ", rdstate: " << ascii_in->rdstate() << endl;
          }
          fileclose();
        }
        else
        {
          if (Verbosity() > 0)
          {
	    cout << "Fun4AllHepMCPileupInputManager::run::" << Name();
	    if (skip) cout << " skip";
            cout << " hepmc evt no: " << evt->event_number() << endl;
          }
          events_total++;
          events_thisfile++;
          // check if the local SubsysReco discards this event
          if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
          {
            ResetEvent();
            //	goto readagain;
          }
          else
          {
            break;  // got the evt, move on
          }
        }
      }  // loop until retrieve a valid event
      if (!skip)
      {
      PHHepMCGenEventMap *geneventmap = hepmc_helper.get_geneventmap();
      PHHepMCGenEvent *genevent = nullptr;
      if (hepmc_helper.get_embedding_id() > 0)
      {
        //! If set_embedding_id(i) with a positive number, the pile up event will be inserted with increasing positive embedding_id. This would be a strange way to use pile up.

        genevent = geneventmap->insert_active_event();
      }
      else
      {
        //! If set_embedding_id(i) with a negative number or 0, the pile up event will be inserted with increasing positive embedding_id. This is the default operation mode.

        genevent = geneventmap->insert_background_event();
      }
      assert(genevent);
      assert(evt);
      genevent->addEvent(evt);
      cout << "handling event " << evt->event_number() << endl;
      hepmc_helper.move_vertex(genevent);
      // place to the crossing center in time
      genevent->moveVertex(0, 0, 0, t0);
      }
    }  //    for (int icollision = 0; icollision < ncollisions; ++icollision)

  }  //  for (int icrossing = _min_crossing; icrossing <= _max_crossing; ++icrossing)

  return 0;
}

int Fun4AllHepMCPileupInputManager::ResetEvent()
{
  m_SignalEventNumber = 0;
  return 0;
}
