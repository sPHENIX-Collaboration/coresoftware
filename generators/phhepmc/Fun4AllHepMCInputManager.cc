#include "Fun4AllHepMCInputManager.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <PHHepMCGenEvent.h>
#include <PHHepMCGenEventMap.h>

#include <ffaobjects/RunHeader.h>

#include <frog/FROG.h>

#include <phool/PHRandomSeed.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>

#include <HepMC/GenEvent.h>
#include <HepMC/IO_GenEvent.h>

#include <TPRegexp.h>
#include <TString.h>


#include <gsl/gsl_const.h>

#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <fstream>
#include <iostream>
#include <istream>
#include <memory>
#include <sstream>
#include <cstdlib>

using namespace std;

static const double toMM = 1.e-12;

Fun4AllHepMCInputManager::Fun4AllHepMCInputManager(const string &name, const string &nodename, const string &topnodename)
  : Fun4AllInputManager(name, nodename, topnodename)
  , isopen(0)
  , events_total(0)
  , events_thisfile(0)
  , readoscar(0)
  , topNodeName(topnodename)
  , ascii_in(nullptr)
  , evt(nullptr)
  , save_evt(nullptr)
  , filestream(nullptr)
  , unzipstream(nullptr)
{
  set_embedding_id(0);  // default embedding ID. Welcome to change via macro

  Fun4AllServer *se = Fun4AllServer::instance();
  topNode = se->topNode(topNodeName.c_str());
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = se->getNode(InputNode(), topNodeName);

  PHHepMCGenEventMap *geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    geneventmap = new PHHepMCGenEventMap();
    PHIODataNode<PHObject> *newmapnode = new PHIODataNode<PHObject>(geneventmap, "PHHepMCGenEventMap", "PHObject");
    dstNode->addNode(newmapnode);
  }

  hepmc_helper.set_geneventmap(geneventmap);

  return;
}

Fun4AllHepMCInputManager::~Fun4AllHepMCInputManager()
{
  fileclose();

  delete ascii_in;
  delete filestream;
  delete unzipstream;
}

int Fun4AllHepMCInputManager::fileopen(const string &filenam)
{
  if (!mySyncManager)
  {
    cout << "Call fileopen only after you registered your Input Manager " << Name() << " with the Fun4AllServer" << endl;
    exit(1);
  }
  if (isopen)
  {
    cout << "Closing currently open file "
         << filename
         << " and opening " << filenam << endl;
    fileclose();
  }
  filename = filenam;
  FROG frog;
  string fname(frog.location(filename.c_str()));
  if (Verbosity() > 0)
  {
    cout << Name() << ": opening file " << fname << endl;
  }

  if (readoscar)
  {
    theOscarFile.open(fname.c_str());
  }
  else
  {
    TString tstr(fname);
    TPRegexp bzip_ext(".bz2$");
    TPRegexp gzip_ext(".gz$");
    if (tstr.Contains(bzip_ext))
    {
      // use boost iosteam library to decompress bz2 on the fly
      filestream = new ifstream(fname.c_str(), std::ios::in | std::ios::binary);
      zinbuffer.push(boost::iostreams::bzip2_decompressor());
      zinbuffer.push(*filestream);
      unzipstream = new istream(&zinbuffer);
      ascii_in = new HepMC::IO_GenEvent(*unzipstream);
    }
    else if (tstr.Contains(gzip_ext))
    {
      // use boost iosream to decompress the gzip file on the fly
      filestream = new ifstream(fname.c_str(), std::ios::in | std::ios::binary);
      zinbuffer.push(boost::iostreams::gzip_decompressor());
      zinbuffer.push(*filestream);
      unzipstream = new istream(&zinbuffer);
      ascii_in = new HepMC::IO_GenEvent(*unzipstream);
    }
    else
    {
      // expects normal ascii hepmc file
      ascii_in = new HepMC::IO_GenEvent(fname, std::ios::in);
    }
  }

  recoConsts *rc = recoConsts::instance();
  static bool run_number_forced = rc->FlagExist("RUNNUMBER");
  if (run_number_forced)
  {
    mySyncManager->CurrentRun(rc->get_IntFlag("RUNNUMBER"));
  }
  else
  {
    mySyncManager->CurrentRun(-1);
  }
  events_thisfile = 0;
  isopen = 1;
  AddToFileOpened(fname);  // add file to the list of files which were opened
  return 0;
}

int Fun4AllHepMCInputManager::run(const int nevents)
{
  // attempt to retrieve a valid event from inputs
  while (true)
  {
    if (!isopen)
    {
      if (FileListEmpty())

      {
        if (Verbosity() > 0)
        {
          cout << "Fun4AllHepMCInputManager::run::" << Name() << ": No Input file open" << endl;
        }
        return -1;
      }
      else
      {
        if (OpenNextFile())
        {
          cout << "Fun4AllHepMCInputManager::run::" << Name() << ": No Input file from filelist opened" << endl;
          return -1;
        }
      }
    }

    if (save_evt)  // if an event was pushed back, copy saved pointer and reset save_evt pointer
    {
      evt = save_evt;
      save_evt = nullptr;
    }
    else
    {
      if (readoscar)
      {
        evt = ConvertFromOscar();
      }
      else
      {
        evt = ascii_in->read_next_event();
      }
    }

    if (!evt)
    {
      if (Verbosity() > 1)
      {
        cout << "Fun4AllHepMCInputManager::run::" << Name()
             << ": error type: " << ascii_in->error_type()
             << ", rdstate: " << ascii_in->rdstate() << endl;
      }
      fileclose();
    }
    else
    {
      mySyncManager->CurrentEvent(evt->event_number());
      if (Verbosity() > 0)
      {
        cout << "Fun4AllHepMCInputManager::run::" << Name()
             << ": hepmc evt no: " << evt->event_number() << endl;
      }

      PHHepMCGenEventMap::Iter ievt =
          hepmc_helper.get_geneventmap()->find(hepmc_helper.get_embedding_id());
      if (ievt != hepmc_helper.get_geneventmap()->end())
      {
        // override existing event
        ievt->second->addEvent(evt);
      }
      else
        hepmc_helper.insert_event(evt);

      events_total++;
      events_thisfile++;

      // check if the local SubsysReco discards this event
      if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
      {
        ResetEvent();
      }
      else
        break;  // have a good event, move on
    }
  }  // attempt to retrieve a valid event from inputs

  return 0;
}

int Fun4AllHepMCInputManager::fileclose()
{
  if (!isopen)
  {
    cout << Name() << ": fileclose: No Input file open" << endl;
    return -1;
  }
  if (readoscar)
  {
    theOscarFile.close();
  }
  else
  {
    delete ascii_in;
    ascii_in = nullptr;
  }
  isopen = 0;
  // if we have a file list, move next entry to top of the list
  // or repeat the same entry again
  UpdateFileList();
  return 0;
}

void Fun4AllHepMCInputManager::Print(const string &what) const
{
  Fun4AllInputManager::Print(what);
  return;
}

int Fun4AllHepMCInputManager::PushBackEvents(const int i)
{
  // PushBackEvents is supposedly pushing events back on the stack which works
  // easily with root trees (just grab a different entry) but hard in these HepMC ASCII files.
  // A special case is when the synchronization fails and we need to only push back a single
  // event. In this case we save the evt pointer as save_evt which is used in the run method
  // instead of getting the next event.
  if (i > 0)
  {
    if (i == 1 && evt)  // check on evt pointer makes sure it is not done from the cmd line
    {
      if (Verbosity() > 3)
      {
        cout << Name() << ": pushing back evt no " << evt->event_number() << endl;
      }
      save_evt = evt;
      return 0;
    }
    cout << PHWHERE << Name()
         << " Fun4AllHepMCInputManager cannot pop back events into file"
         << endl;
    return -1;
  }
  if (!isopen)
  {
    cout << PHWHERE << Name()
         << " no file opened yet" << endl;
    return -1;
  }
  // Skipping events is implemented as
  // pushing a negative number of events on the stack, so in order to implement
  // the skipping of events we read -i events.
  int nevents = -i;  // negative number of events to push back -> skip num events
  int errorflag = 0;
  while (nevents > 0 && !errorflag)
  {
    evt = ascii_in->read_next_event();
    if (!evt)
    {
      cout << "Error after skipping " << i - nevents << endl;
      cout << "error type: " << ascii_in->error_type()
           << ", rdstate: " << ascii_in->rdstate() << endl;
      errorflag = -1;
      fileclose();
    }
    else
    {
      if (Verbosity() > 3)
      {
        cout << "Skipping evt no: " << evt->event_number() << endl;
      }
    }
    delete evt;
    nevents--;
  }
  return errorflag;
}

HepMC::GenEvent *
Fun4AllHepMCInputManager::ConvertFromOscar()
{
  delete evt;
  evt = nullptr;
  if (theOscarFile.eof())  // if the file is exhausted bail out during this next read
  {
    return evt;
  }

  //use PHENIX unit
  evt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::CM);

  if (Verbosity() > 1) cout << "Reading Oscar Event " << events_total << endl;
  //Grab New Event From Oscar
  string theLine;
  vector<vector<double> > theEventVec;
  vector<HepMC::FourVector> theVtxVec;
  while (getline(theOscarFile, theLine))
  {
    if (theLine.find("#") == 0) continue;
    vector<double> theInfo;  //format: N,pid,px,py,pz,E,mass,xvtx,yvtx,zvtx,?
    double number;
    for (istringstream numbers_iss(theLine); numbers_iss >> number;)
    {
      theInfo.push_back(number);
    }

    if (theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] == 0)
    {
      break;
    }
    else if (theInfo.size() == 2 && theInfo[0] == 0 && theInfo[1] > 0)
    {
      continue;
    }
    else
    {
      theEventVec.push_back(theInfo);
      HepMC::FourVector vert(theInfo[8] * toMM, theInfo[9] * toMM, theInfo[10] * toMM, theInfo[11]);
      theVtxVec.push_back(vert);
    }

  }  //while(getline)

  //Set Event Number
  evt->set_event_number(events_total);

  //Loop Over One Event, Fill HepMC
  for (unsigned int i = 0; i < theEventVec.size(); i++)
  {
    //int N = (int)theEventVec[i][0];
    int pid = (int) theEventVec[i][1];
    double px = theEventVec[i][3];
    double py = theEventVec[i][4];
    double pz = theEventVec[i][5];
    double E = theEventVec[i][6];
    double m = theEventVec[i][7];
    int status = 1;  //oscar only writes final state particles

    HepMC::GenVertex *v = new HepMC::GenVertex(theVtxVec[i]);
    evt->add_vertex(v);

    HepMC::GenParticle *p = new HepMC::GenParticle(HepMC::FourVector(px, py, pz, E), pid, status);
    p->setGeneratedMass(m);
    p->suggest_barcode(i + 1);
    v->add_particle_out(p);
  }
  if (Verbosity() > 3)
  {
    evt->print();
  }
  return evt;
}
