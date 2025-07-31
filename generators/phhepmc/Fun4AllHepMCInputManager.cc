#include "Fun4AllHepMCInputManager.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"

#include <frog/FROG.h>

#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInpu...
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>  // for GenParticle
#include <HepMC/GenVertex.h>    // for GenVertex
#include <HepMC/IO_GenEvent.h>
#include <HepMC/SimpleVector.h>  // for FourVector
#include <HepMC/Units.h>         // for CM, GEV

#include <TPRegexp.h>
#include <TString.h>

#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>  // for _Rb_tree_it...
#include <sstream>
#include <vector>  // for vector

static const double toMM = 1.e-12;

Fun4AllHepMCInputManager::Fun4AllHepMCInputManager(const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllInputManager(name, nodename, topnodename)
  , topNodeName(topnodename)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  topNode = se->topNode(topNodeName);
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = se->getNode(InputNode(), topNodeName);

  PHHepMCGenEventMap *geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    geneventmap = new PHHepMCGenEventMap();
    PHIODataNode<PHObject> *newmapnode = new PHIODataNode<PHObject>(geneventmap, "PHHepMCGenEventMap", "PHObject");
    dstNode->addNode(newmapnode);
  }

  PHHepMCGenHelper::set_geneventmap(geneventmap);

  return;
}

Fun4AllHepMCInputManager::~Fun4AllHepMCInputManager()
{
  fileclose();
  if (!m_HepMCTmpFile.empty())
  {
    // okay if the file does not exist
    remove(m_HepMCTmpFile.c_str());
  }
  delete ascii_in;
  delete filestream;
  delete unzipstream;
}

int Fun4AllHepMCInputManager::fileopen(const std::string &filenam)
{
  if (!MySyncManager())
  {
    std::cout << "Call fileopen only after you registered your Input Manager " << Name() << " with the Fun4AllServer" << std::endl;
    exit(1);
  }
  if (IsOpen())
  {
    std::cout << "Closing currently open file "
              << filename
              << " and opening " << filenam << std::endl;
    fileclose();
  }
  filename = filenam;
  FROG frog;
  std::string fname(frog.location(filename));
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << fname << std::endl;
  }

  if (m_ReadOscarFlag)
  {
    theOscarFile.open(fname);
  }
  else
  {
    TString tstr(fname);
    TPRegexp bzip_ext(".bz2$");
    TPRegexp gzip_ext(".gz$");
    if (tstr.Contains(bzip_ext))
    {
      // use boost iosteam library to decompress bz2 on the fly
      filestream = new std::ifstream(fname, std::ios::in | std::ios::binary);
      zinbuffer.push(boost::iostreams::bzip2_decompressor());
      zinbuffer.push(*filestream);
      unzipstream = new std::istream(&zinbuffer);
      ascii_in = new HepMC::IO_GenEvent(*unzipstream);
    }
    else if (tstr.Contains(gzip_ext))
    {
      // use boost iosream to decompress the gzip file on the fly
      filestream = new std::ifstream(fname, std::ios::in | std::ios::binary);
      zinbuffer.push(boost::iostreams::gzip_decompressor());
      zinbuffer.push(*filestream);
      unzipstream = new std::istream(&zinbuffer);
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
    MySyncManager()->CurrentRun(rc->get_IntFlag("RUNNUMBER"));
  }
  else
  {
    MySyncManager()->CurrentRun(-1);
  }
  events_thisfile = 0;
  IsOpen(1);
  AddToFileOpened(fname);  // add file to the list of files which were opened
  return 0;
}

int Fun4AllHepMCInputManager::run(const int /*nevents*/)
{
  // attempt to retrieve a valid event from inputs
  while (true)
  {
    if (!IsOpen())
    {
      if (FileListEmpty())

      {
        if (Verbosity() > 0)
        {
          std::cout << "Fun4AllHepMCInputManager::run::" << Name() << ": No Input file open" << std::endl;
        }
        return -1;
      }
      else
      {
        if (OpenNextFile())
        {
          std::cout << "Fun4AllHepMCInputManager::run::" << Name() << ": No Input file from filelist opened" << std::endl;
          return -1;
        }
      }
    }

    if (m_EventPushedBackFlag)  // if an event was pushed back, reuse save copy
    {
      HepMC::IO_GenEvent ascii_tmp_in(m_HepMCTmpFile, std::ios::in);
      ascii_tmp_in >> evt;
      m_EventPushedBackFlag = 0;
      remove(m_HepMCTmpFile.c_str());
    }
    else
    {
      if (m_ReadOscarFlag)
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
        std::cout << "Fun4AllHepMCInputManager::run::" << Name()
                  << ": error type: " << ascii_in->error_type()
                  << ", rdstate: " << ascii_in->rdstate() << std::endl;
      }
      fileclose();
    }
    else
    {
      MySyncManager()->CurrentEvent(evt->event_number());
      if (Verbosity() > 0)
      {
        std::cout << "Fun4AllHepMCInputManager::run::" << Name()
                  << ": hepmc evt no: " << evt->event_number() << std::endl;
      }
      m_MyEvent.push_back(evt->event_number());
      PHHepMCGenEventMap::Iter ievt = PHHepMCGenHelper::get_geneventmap()->find(PHHepMCGenHelper::get_embedding_id());
      if (ievt != PHHepMCGenHelper::get_geneventmap()->end())
      {
        // override existing event
        ievt->second->addEvent(evt);
      }
      else
      {
        PHHepMCGenHelper::insert_event(evt);
      }

      events_total++;
      events_thisfile++;

      // check if the local SubsysReco discards this event
      if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK)
      {
        // if this event is discarded we only need to remove the event from the list of event numbers
        // the new event will overwrite the event on the node tree without issues
        m_MyEvent.pop_back();
      }
      else
      {
        break;  // have a good event, move on
      }
    }
  }  // attempt to retrieve a valid event from inputs
  return Fun4AllReturnCodes::EVENT_OK;
}

int Fun4AllHepMCInputManager::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  if (m_ReadOscarFlag)
  {
    theOscarFile.close();
  }
  else
  {
    delete ascii_in;
    ascii_in = nullptr;
  }
  IsOpen(0);
  // if we have a file list, move next entry to top of the list
  // or repeat the same entry again
  UpdateFileList();
  return 0;
}

void Fun4AllHepMCInputManager::Print(const std::string &what) const
{
  Fun4AllInputManager::Print(what);
  std::cout << Name() << " Vertex Settings: " << std::endl;
  PHHepMCGenHelper::Print(what);
  return;
}

int Fun4AllHepMCInputManager::PushBackEvents(const int i)
{
  // PushBackEvents is supposedly pushing events back on the stack which works
  // easily with root trees (just grab a different entry) but hard in these HepMC ASCII files.
  // A special case is when the synchronization fails and we need to only push back a single
  // event. In this case we save the evt in a temporary file which is read back in the run method
  // instead of getting the next event.
  if (i > 0)
  {
    if (i == 1 && evt)  // check on evt pointer makes sure it is not done from the cmd line
    {
      // root barfs when writing the node to the output.
      // Saving the pointer - even using a deep copy and reusing it did not work
      // The hackaround which works is to write this event into a temporary file and read it back
      if (m_HepMCTmpFile.empty())
      {
        // we need to create this filename just once, we reuse it. Do it only if we need it
        m_HepMCTmpFile = "/tmp/HepMCTmpEvent-" + Name() + "-" + std::to_string(getpid()) + ".hepmc";
      }
      HepMC::IO_GenEvent ascii_io(m_HepMCTmpFile, std::ios::out);
      ascii_io << evt;
      m_EventPushedBackFlag = 1;
      m_MyEvent.pop_back();

      if (Verbosity() > 3)
      {
        std::cout << Name() << ": pushing back evt no " << evt->event_number() << std::endl;
      }
      return 0;
    }
    std::cout << PHWHERE << Name()
              << " Fun4AllHepMCInputManager cannot pop back more than 1 event"
              << std::endl;
    return -1;
  }
  if (!IsOpen())
  {
    std::cout << PHWHERE << Name()
              << " no file opened yet" << std::endl;
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
      std::cout << "Error after skipping " << i - nevents << std::endl;
      std::cout << "error type: " << ascii_in->error_type()
                << ", rdstate: " << ascii_in->rdstate() << std::endl;
      errorflag = -1;
      fileclose();
    }
    else
    {
      m_MyEvent.push_back(evt->event_number());
      if (Verbosity() > 3)
      {
        std::cout << "Skipping evt no: " << evt->event_number() << std::endl;
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
  if (theOscarFile.eof())  // if the file is exhausted bail out during this next read
  {
    return nullptr;
  }

  delete evt;
  //use PHENIX unit
  evt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::CM);

  if (Verbosity() > 1) std::cout << "Reading Oscar Event " << events_total << std::endl;
  //Grab New Event From Oscar
  std::string theLine;
  std::vector<std::vector<double> > theEventVec;
  std::vector<HepMC::FourVector> theVtxVec;
  while (getline(theOscarFile, theLine))
  {
    if (theLine.compare(0, 1, "#") == 0) continue;
    std::vector<double> theInfo;  //format: N,pid,px,py,pz,E,mass,xvtx,yvtx,zvtx,?
    double number = NAN;
    for (std::istringstream numbers_iss(theLine); numbers_iss >> number;)
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

int Fun4AllHepMCInputManager::ResetEvent()
{
  m_MyEvent.clear();
  return 0;
}

int Fun4AllHepMCInputManager::MyCurrentEvent(const unsigned int index) const
{
  if (m_MyEvent.empty())
  {
    return 0;
  }
  return m_MyEvent.at(index);
}

