#include "Fun4AllServer.h"

#include "Fun4AllHistoBinDefs.h"
#include "Fun4AllHistoManager.h"  // for Fun4AllHistoManager
#include "Fun4AllMemoryTracker.h"
#include "Fun4AllMonitoring.h"
#include "Fun4AllOutputManager.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllSyncManager.h"
#include "SubsysReco.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHNodeReset.h>
#include <phool/PHObject.h>
#include <phool/PHPointerListIterator.h>
#include <phool/PHTimeStamp.h>
#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <Rtypes.h>  // for kMAXSIGNALS
#include <TDirectory.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSysEvtHandler.h>  // for ESignals

#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>  // for allocator_traits<>::value_type
#include <sstream>

//#define FFAMEMTRACKER

Fun4AllServer *Fun4AllServer::__instance = nullptr;

Fun4AllServer *Fun4AllServer::instance()
{
  if (__instance)
  {
    return __instance;
  }
  __instance = new Fun4AllServer();
  return __instance;
}

Fun4AllServer::Fun4AllServer(const std::string &name)
  : Fun4AllBase(name)
#ifdef FFAMEMTRACKER
  , ffamemtracker(Fun4AllMemoryTracker::instance())
#endif
{
  InitAll();
  return;
}

Fun4AllServer::~Fun4AllServer()
{
  Reset();
  delete beginruntimestamp;
  while (Subsystems.begin() != Subsystems.end())
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      Subsystems.back().first->Verbosity(Verbosity());
    }
    delete Subsystems.back().first;
    Subsystems.pop_back();
  }
  while (HistoManager.begin() != HistoManager.end())
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      HistoManager.back()->Verbosity(Verbosity());
    }
    delete HistoManager.back();
    HistoManager.pop_back();
  }
  while (OutputManager.begin() != OutputManager.end())
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      OutputManager.back()->Verbosity(Verbosity());
    }
    delete OutputManager.back();
    OutputManager.pop_back();
  }
  while (SyncManagers.begin() != SyncManagers.end())
  {
    SyncManagers.back()->Verbosity(Verbosity());
    delete SyncManagers.back();
    SyncManagers.pop_back();
  }
  while (topnodemap.begin() != topnodemap.end())
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      topnodemap.begin()->second->print();
    }
    delete topnodemap.begin()->second;
    topnodemap.erase(topnodemap.begin());
  }
  while (TDirCollection.begin() != TDirCollection.end())
  {
    delete TDirCollection.back();
    TDirCollection.pop_back();
  }
  recoConsts *rc = recoConsts::instance();
  delete rc;
  delete ffamemtracker;
  __instance = nullptr;
  return;
}

void Fun4AllServer::InitAll()
{
  // first remove stupid root signal handler to get
  // decent crashes with debuggable core file
  for (int i = 0; i < kMAXSIGNALS; i++)
  {
    gSystem->IgnoreSignal((ESignals) i);
  }
  Fun4AllMonitoring::instance()->Snapshot("StartUp");
  std::string histomanagername;
  histomanagername = Name() + "HISTOS";
  ServerHistoManager = new Fun4AllHistoManager(histomanagername);
  registerHistoManager(ServerHistoManager);
  double uplim = NFRAMEWORKBINS - 0.5;
  FrameWorkVars = new TH1D("FrameWorkVars", "FrameWorkVars", NFRAMEWORKBINS, -0.5, uplim);
  registerHisto("FrameWorkVars", FrameWorkVars);
  defaultSyncManager = new Fun4AllSyncManager("DefaultSyncManager");
  SyncManagers.push_back(defaultSyncManager);
  TopNode = new PHCompositeNode("TOP");
  topnodemap["TOP"] = TopNode;
  default_Tdirectory = gDirectory->GetPath();
  InitNodeTree(TopNode);
  return;
}

int Fun4AllServer::dumpHistos(const std::string &filename, const std::string &openmode)
{
  int iret = 0;
  std::cout << "Fun4AllServer::dumpHistos() dumping histograms" << std::endl;
  if (!filename.empty())
  {
    ServerHistoManager->setOutfileName(filename);
  }
  std::vector<Fun4AllHistoManager *>::const_iterator hiter;
  for (hiter = HistoManager.begin(); hiter != HistoManager.end(); ++hiter)
  {
    iret += (*hiter)->dumpHistos("", openmode);
  }
  return iret;
}

bool Fun4AllServer::registerHisto(TNamed *h1d, const int replace)
{
  return ServerHistoManager->registerHisto(h1d, replace);
}

bool Fun4AllServer::registerHisto(const std::string &hname, TNamed *h1d, const int replace)
{
  return ServerHistoManager->registerHisto(hname, h1d, replace);
}

int Fun4AllServer::isHistoRegistered(const std::string &name) const
{
  int iret = ServerHistoManager->isHistoRegistered(name);
  return iret;
}

int Fun4AllServer::registerSubsystem(SubsysReco *subsystem, const std::string &topnodename)
{
  Fun4AllServer *se = Fun4AllServer::instance();

  // if somebody opens a TFile (or changes the gDirectory) in the ctor
  // we need to set it to a "known" directory
  gROOT->cd(default_Tdirectory.c_str());
  std::string currdir = gDirectory->GetPath();
  TDirectory *tmpdir = gDirectory;
  if (!tmpdir->FindObject(topnodename.c_str()))
  {
    tmpdir = tmpdir->mkdir(topnodename.c_str());
    if (!tmpdir)
    {
      std::cout << "Error creating TDirectory topdir " << topnodename.c_str() << std::endl;
      exit(1);
    }
    // store the TDir pointer so it can be cleaned up in the dtor
    // if one deletes it here the Histograms are dangling somewhere
    // in root (at least according to valgrind the delete doesn't work
    // properly anymore)
    TDirCollection.push_back(tmpdir);
  }
  gROOT->cd(topnodename.c_str());
  tmpdir = gDirectory;
  if (!tmpdir->FindObject(subsystem->Name().c_str()))
  {
    tmpdir = tmpdir->mkdir(subsystem->Name().c_str());
    if (!tmpdir)
    {
      std::cout << "Error creating TDirectory subdir " << subsystem->Name() << std::endl;
      exit(1);
    }
    // store the TDir pointer so it can be cleaned up in the dtor
    // if one deletes it here the Histograms are dangling somewhere
    // in root
    TDirCollection.push_back(tmpdir);
  }
  PHCompositeNode *subsystopNode = se->topNode(topnodename);
  std::pair<SubsysReco *, PHCompositeNode *> newsubsyspair(subsystem, subsystopNode);
  int iret = 0;
  try
  {
    std::string memory_tracker_name = subsystem->Name() + "_" + topnodename;
#ifdef FFAMEMTRACKER
    ffamemtracker->Start(memory_tracker_name, "SubsysReco");
#endif
    iret = subsystem->Init(subsystopNode);
#ifdef FFAMEMTRACKER
    ffamemtracker->Stop(memory_tracker_name, "SubsysReco");
#endif
  }
  catch (const std::exception &e)
  {
    std::cout << PHWHERE << " caught exception thrown during SubsysReco::Init() from "
              << subsystem->Name() << std::endl;
    std::cout << "error: " << e.what() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::Init() from "
              << subsystem->Name() << std::endl;
    exit(1);
  }
  gROOT->cd(currdir.c_str());
  if (iret)
  {
    if (iret == Fun4AllReturnCodes::DONOTREGISTERSUBSYSTEM)
    {
      if (Verbosity() >= VERBOSITY_SOME)
      {
        std::cout << "Not Registering Subsystem " << subsystem->Name() << std::endl;
      }
      return 0;
    }
    std::cout << PHWHERE << " Error from Init() method by "
              << subsystem->Name() << ", return code: " << iret << std::endl;
    return iret;
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Registering Subsystem " << subsystem->Name() << std::endl;
  }
  Subsystems.push_back(newsubsyspair);
  std::string timer_name;
  timer_name = subsystem->Name() + "_" + topnodename;
  PHTimer timer(timer_name);
  if (timer_map.find(timer_name) == timer_map.end())
  {
    timer_map.insert(make_pair(timer_name, timer));
  }
  RetCodes.push_back(iret);  // vector with return codes
  return 0;
}

int Fun4AllServer::unregisterSubsystem(SubsysReco *subsystem)
{
  std::pair<SubsysReco *, PHCompositeNode *> subsyspair(subsystem, 0);
  DeleteSubsystems.push_back(subsyspair);
  unregistersubsystem = 1;
  return 0;
}

int Fun4AllServer::unregisterSubsystemsNow()
{
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::iterator sysiter, removeiter;
  for (removeiter = DeleteSubsystems.begin();
       removeiter != DeleteSubsystems.end();
       ++removeiter)
  {
    int index = 0;
    int foundit = 0;
    for (sysiter = Subsystems.begin(); sysiter != Subsystems.end(); ++sysiter)
    {
      if ((*sysiter).first == (*removeiter).first)
      {
        foundit = 1;
        break;
      }
      index++;
    }
    if (!foundit)
    {
      std::cout << "unregisterSubsystem: Could not find SubsysReco "
                << (*removeiter).first->Name()
                << " in Fun4All Reco Module list" << std::endl;
      delete (*removeiter).first;
      continue;
    }
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Removing Subsystem: " << (*removeiter).first->Name()
                << " at index " << index << std::endl;
    }
    Subsystems.erase(Subsystems.begin() + index);
    delete (*removeiter).first;
    // also update the vector with return codes
    RetCodes.erase(RetCodes.begin() + index);
    std::vector<Fun4AllOutputManager *>::iterator outiter;
    for (outiter = OutputManager.begin(); outiter != OutputManager.end(); ++outiter)
    {
      UpdateEventSelector(*outiter);
    }
  }
  unregistersubsystem = 0;
  DeleteSubsystems.clear();
  return 0;
}

SubsysReco *
Fun4AllServer::getSubsysReco(const std::string &name)
{
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::iterator sysiter;
  for (sysiter = Subsystems.begin(); sysiter != Subsystems.end(); ++sysiter)
  {
    if ((*sysiter).first->Name() == name)
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "Found Subsystem " << name << std::endl;
      }
      return (*sysiter).first;
    }
  }
  std::cout << "Could not find SubsysReco " << name << std::endl;
  return nullptr;
}

int Fun4AllServer::AddComplaint(const std::string &complaint, const std::string &remedy)
{
  ScreamEveryEvent++;
  std::string separatorstring = "------------------------------";
  std::ostringstream complaintno;
  complaintno << "Problem No " << ScreamEveryEvent;

  ComplaintList.push_back(separatorstring);
  ComplaintList.push_back(complaintno.str());
  ComplaintList.push_back(complaint);
  ComplaintList.emplace_back(" ");
  ComplaintList.emplace_back("Remedy:");
  ComplaintList.push_back(remedy);
  ComplaintList.push_back(separatorstring);
  return 0;
}

int Fun4AllServer::registerOutputManager(Fun4AllOutputManager *manager)
{
  std::vector<Fun4AllOutputManager *>::iterator iter;
  for (iter = OutputManager.begin(); iter != OutputManager.end(); ++iter)
  {
    if ((*iter)->Name() == manager->Name())
    {
      std::cout << "OutputManager " << manager->Name() << " allready in list" << std::endl;
      return -1;
    }
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Registering OutputManager " << manager->Name() << std::endl;
  }
  UpdateEventSelector(manager);
  OutputManager.push_back(manager);
  return 0;
}

int Fun4AllServer::UpdateEventSelector(Fun4AllOutputManager *manager)
{
  std::vector<std::string>::iterator striter;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::const_iterator subsysiter;

tryagain:
  manager->RecoModuleIndex()->clear();
  for (striter = manager->EventSelector()->begin(); striter != manager->EventSelector()->end(); ++striter)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      std::cout << PHWHERE << "striter: " << *striter << std::endl;
    }
    unsigned index = 0;
    int found = 0;
    for (subsysiter = Subsystems.begin(); subsysiter != Subsystems.end(); ++subsysiter)
    {
      if (*striter == (*subsysiter).first->Name())
      {
        manager->RecoModuleIndex()->push_back(index);
        if (Verbosity() >= VERBOSITY_EVEN_MORE)
        {
          std::cout << PHWHERE << "setting RecoModuleIndex to " << index << std::endl;
        }
        found = 1;
        break;
      }
      index++;
    }
    if (!found)
    {
      std::cout << "Could not find module " << *striter
                << ", removing it from list of event selector modules" << std::endl;
      manager->EventSelector()->erase(striter);
      goto tryagain;
    }
  }
  return 0;
}

Fun4AllOutputManager *
Fun4AllServer::getOutputManager(const std::string &name)
{
  std::vector<Fun4AllOutputManager *>::iterator iter;
  for (iter = OutputManager.begin(); iter != OutputManager.end(); ++iter)
  {
    if (name == (*iter)->Name())
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "Found OutputManager " << name << std::endl;
      }
      return *iter;
    }
  }
  std::cout << "Could not find OutputManager" << name << std::endl;
  return nullptr;
}

Fun4AllHistoManager *
Fun4AllServer::getHistoManager(const std::string &name)
{
  std::vector<Fun4AllHistoManager *>::iterator iter;
  for (iter = HistoManager.begin(); iter != HistoManager.end(); ++iter)
  {
    if ((*iter)->Name() == name)
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "Found HistoManager " << name << std::endl;
      }
      return *iter;
    }
  }
  if (Verbosity() >= VERBOSITY_MORE)
  {
    std::cout << "Could not find HistoManager " << name << std::endl;
  }
  return nullptr;
}

int Fun4AllServer::registerHistoManager(Fun4AllHistoManager *manager)
{
  std::vector<Fun4AllHistoManager *>::iterator iter;
  for (iter = HistoManager.begin(); iter != HistoManager.end(); ++iter)
  {
    if ((*iter)->Name() == manager->Name())
    {
      std::cout << "HistoManager " << manager->Name() << " allready in list" << std::endl;
      return -1;
    }
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Registering HistoManager " << manager->Name() << std::endl;
  }
  HistoManager.push_back(manager);
  return 0;
}

TNamed *
Fun4AllServer::getHisto(const unsigned int ihisto) const
{
  return ServerHistoManager->getHisto(ihisto);
}

std::string
Fun4AllServer::getHistoName(const unsigned int ihisto) const
{
  return (ServerHistoManager->getHistoName(ihisto));
}

TNamed *Fun4AllServer::getHisto(const std::string &hname) const
{
  return (ServerHistoManager->getHisto(hname));
}

int Fun4AllServer::process_event()
{
  eventcounter++;
  unsigned icnt = 0;
  int eventbad = 0;
  if (ScreamEveryEvent)
  {
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "Now that I have your attention, please fix the following "
              << ScreamEveryEvent << " problem(s):" << std::endl;
    std::vector<std::string>::const_iterator viter;
    for (viter = ComplaintList.begin(); viter != ComplaintList.end(); ++viter)
    {
      std::cout << *viter << std::endl;
    }
    std::cout << " " << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
  }
  if (unregistersubsystem)
  {
    unregisterSubsystemsNow();
  }
  gROOT->cd(default_Tdirectory.c_str());
  std::string currdir = gDirectory->GetPath();
  for (auto & Subsystem : Subsystems)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Fun4AllServer::process_event processing " << Subsystem.first->Name() << std::endl;
    }
    std::ostringstream newdirname;
    newdirname << Subsystem.second->getName() << "/" << Subsystem.first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
                << Subsystem.second->getName()
                << " - send e-mail to off-l with your macro" << std::endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "process_event: cded to " << newdirname.str().c_str() << std::endl;
      }
    }

    try
    {
      std::string timer_name;
      timer_name = Subsystem.first->Name() + "_" + Subsystem.second->getName();
      std::map<const std::string, PHTimer>::iterator titer = timer_map.find(timer_name);
      bool timer_found = false;
      if (titer != timer_map.end())
      {
        timer_found = true;
        titer->second.restart();
      }
      else
      {
        std::cout << "could not find timer for " << timer_name << std::endl;
      }
#ifdef FFAMEMTRACKER
      ffamemtracker->Start(timer_name, "SubsysReco");
      ffamemtracker->Snapshot("Fun4AllServerProcessEvent");
#endif
      int retcode = Subsystem.first->process_event(Subsystem.second);
#ifdef FFAMEMTRACKER
      ffamemtracker->Snapshot("Fun4AllServerProcessEvent");
#endif
      // we have observed an index overflow in RetCodes. I assume it is some
      // memory corruption elsewhere which hits the icnt variable. Rather than
      // the previous [], use at() which does bounds checking and throws an
      // exception which will allow us to catch this and print out icnt and the size
      try
      {
        RetCodes.at(icnt) = retcode;
      }
      catch (const std::exception &e)
      {
        std::cout << PHWHERE << " caught exception thrown during RetCodes.at(icnt)" << std::endl;
        std::cout << "RetCodes.size(): " << RetCodes.size() << ", icnt: " << icnt << std::endl;
        std::cout << "error: " << e.what() << std::endl;
        gSystem->Exit(1);
      }
      if (timer_found)
      {
        titer->second.stop();
      }
#ifdef FFAMEMTRACKER
      ffamemtracker->Stop(timer_name, "SubsysReco");
#endif
    }
    catch (const std::exception &e)
    {
      std::cout << PHWHERE << " caught exception thrown during process_event from "
                << Subsystem.first->Name() << std::endl;
      std::cout << "error: " << e.what() << std::endl;
      gSystem->Exit(1);
    }
    catch (...)
    {
      std::cout << PHWHERE << " caught unknown type exception thrown during process_event from "
                << Subsystem.first->Name() << std::endl;
      exit(1);
    }
    if (RetCodes[icnt])
    {
      if (RetCodes[icnt] == Fun4AllReturnCodes::DISCARDEVENT)
      {
        if (Verbosity() >= VERBOSITY_EVEN_MORE)
        {
          std::cout << "Fun4AllServer::Discard Event by " << Subsystem.first->Name() << std::endl;
        }
      }
      else if (RetCodes[icnt] == Fun4AllReturnCodes::ABORTEVENT)
      {
        retcodesmap[Fun4AllReturnCodes::ABORTEVENT]++;
        eventbad = 1;
        if (Verbosity() >= VERBOSITY_MORE)
        {
          std::cout << "Fun4AllServer::Abort Event by " << Subsystem.first->Name() << std::endl;
        }
        break;
      }
      else if (RetCodes[icnt] == Fun4AllReturnCodes::ABORTRUN)
      {
        retcodesmap[Fun4AllReturnCodes::ABORTRUN]++;
        std::cout << "Fun4AllServer::Abort Run by " << Subsystem.first->Name() << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      else
      {
        std::cout << "Fun4AllServer::Unknown return code: "
                  << RetCodes[icnt] << " from process_event method of "
                  << Subsystem.first->Name() << std::endl;
        std::cout << "This smells like an uninitialized return code and" << std::endl;
        std::cout << "it is too dangerous to continue, this Run will be aborted" << std::endl;
        std::cout << "If you do not know how to fix this please send mail to" << std::endl;
        std::cout << "phenix-off-l with this message" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
    icnt++;
  }
  if (!eventbad)
  {
    retcodesmap[Fun4AllReturnCodes::EVENT_OK]++;
  }

  gROOT->cd(currdir.c_str());

  //  mainIter.print();
  if (!OutputManager.empty() && !eventbad)  // there are registered IO managers and
  // the event is not flagged bad
  {
    PHNodeIterator iter(TopNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    if (dstNode)
    {
      // check if we have same number of nodes. After first event is
      // written out root I/O doesn't permit adding nodes, otherwise
      // events get out of sync
      static int first = 1;
      int newcount = CountOutNodes(dstNode);
      if (first)
      {
        first = 0;
        OutNodeCount = newcount;      // save number of nodes before first write
        MakeNodesTransient(dstNode);  // make all nodes transient before 1st write in case someone sneaked a node in at the first event
      }

      if (OutNodeCount != newcount)
      {
        iter.print();
        std::cout << PHWHERE << " FATAL: Someone changed the number of Output Nodes on the fly, from " << OutNodeCount << " to " << newcount << std::endl;
        exit(1);
      }
      std::vector<Fun4AllOutputManager *>::iterator iterOutMan;
      for (iterOutMan = OutputManager.begin(); iterOutMan != OutputManager.end(); ++iterOutMan)
      {
        if (!(*iterOutMan)->DoNotWriteEvent(&RetCodes))
        {
          if (Verbosity() >= VERBOSITY_MORE)
          {
            std::cout << "Writing Event for " << (*iterOutMan)->Name() << std::endl;
          }
#ifdef FFAMEMTRACKER
          ffamemtracker->Snapshot("Fun4AllServerOutputManager");
          ffamemtracker->Start((*iterOutMan)->Name(), "OutputManager");
#endif
          (*iterOutMan)->WriteGeneric(dstNode);
#ifdef FFAMEMTRACKER
          ffamemtracker->Stop((*iterOutMan)->Name(), "OutputManager");
          ffamemtracker->Snapshot("Fun4AllServerOutputManager");
#endif
        }
        else
        {
          if (Verbosity() >= VERBOSITY_MORE)
          {
            std::cout << "Not Writing Event for " << (*iterOutMan)->Name() << std::endl;
          }
        }
      }
    }
  }
  for (auto & Subsystem : Subsystems)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      std::cout << "Fun4AllServer::process_event Resetting Event " << Subsystem.first->Name() << std::endl;
    }
    Subsystem.first->ResetEvent(Subsystem.second);
  }
  for (auto &syncman : SyncManagers)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      std::cout << "Fun4AllServer::process_event Resetting Event for Sync Manager " << syncman->Name() << std::endl;
    }
    syncman->ResetEvent();
  }
  Fun4AllMonitoring::instance()->Snapshot("Event");
  ResetNodeTree();
  return 0;
}

int Fun4AllServer::ResetNodeTree()
{
  std::vector<std::string> ResetNodeList;
  ResetNodeList.emplace_back("DST");
  PHNodeReset reset;
  reset.Verbosity(Verbosity() > 2 ? Verbosity() - 2 : 0);  // one lower verbosity level than Fun4AllServer
  std::map<std::string, PHCompositeNode *>::const_iterator iter;
  for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
  {
    PHNodeIterator mainIter((*iter).second);
    for (const auto & nodename : ResetNodeList)
    {
      if (mainIter.cd(nodename))
      {
        mainIter.forEach(reset);
        mainIter.cd();
      }
    }
  }
  return 0;  // anything except 0 would abort the event loop in pmonitor
}

int Fun4AllServer::Reset()
{
  int i = 0;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::iterator iter;
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      std::cout << "Fun4AllServer::Reset Resetting " << (*iter).first->Name() << std::endl;
    }
    i += (*iter).first->Reset((*iter).second);
  }
  std::vector<Fun4AllHistoManager *>::iterator hiter;
  for (hiter = HistoManager.begin(); hiter != HistoManager.end(); ++hiter)
  {
    (*hiter)->Reset();
  }
  return i;
}

int Fun4AllServer::BeginRunTimeStamp(PHTimeStamp &TimeStp)
{
  beginruntimestamp = new PHTimeStamp(TimeStp);
  std::cout << "Setting BOR timestamp to ";
  beginruntimestamp->print();
  std::cout << std::endl;
  bortime_override = 1;
  return 0;
}

int Fun4AllServer::BeginRun(const int runno)
{
  eventcounter = 0;  // reset event counter for every new run
#ifdef FFAMEMTRACKER
  ffamemtracker->Snapshot("Fun4AllServerBeginRun");
#endif
  if (!bortime_override)
  {
    if (beginruntimestamp)
    {
      delete beginruntimestamp;
    }
    beginruntimestamp = new PHTimeStamp();
  }
  else
  {
    std::cout << "overriding BOR timestamp by ";
    beginruntimestamp->print();
    std::cout << std::endl;
    //rc->set_TimeStamp(*beginruntimestamp);
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Fun4AllServer::BeginRun: Run number " << runno << " uses RECO TIMESTAMP: ";
    beginruntimestamp->print();
    std::cout << std::endl;
  }
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::iterator iter;
  int iret = 0;

  // check if any registered SubsysReco wants to be dropped and
  // remove it from the list before its BeginRun is executed
  if (unregistersubsystem)
  {
    unregisterSubsystemsNow();
  }

  // we have to do the same TDirectory games as in the Init methods
  // save the current dir, cd to the subsystem name dir (which was
  // created in init) call the InitRun of the module and cd back

  gROOT->cd(default_Tdirectory.c_str());
  std::string currdir = gDirectory->GetPath();

  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    iret = BeginRunSubsystem(*iter);
  }
  for (; !NewSubsystems.empty(); NewSubsystems.pop_front())
  {
    registerSubsystem((NewSubsystems.front()).first, (NewSubsystems.front()).second);
    BeginRunSubsystem(std::make_pair(NewSubsystems.front().first, topNode(NewSubsystems.front().second)));
  }
  gROOT->cd(currdir.c_str());
  // disconnect from DB to save resources on DB machine
  // PdbCal leaves the DB connection open (PdbCal will reconnect without
  // problem if neccessary)
  if (!keep_db_connected)
  {
    DisconnectDB();
  }
  else
  {
    std::cout << "WARNING WARNING, DBs will not be disconnected" << std::endl;
    std::cout << "This is for DB server testing purposes only" << std::endl;
    std::cout << "If you do not test our DB servers, remove" << std::endl;
    std::cout << "Fun4AllServer->KeepDBConnection()" << std::endl;
    std::cout << "from your macro" << std::endl;
  }
  // print out all node trees
  Print("NODETREE");
#ifdef FFAMEMTRACKER
  ffamemtracker->Snapshot("Fun4AllServerBeginRun");
#endif
  return iret;
}

int Fun4AllServer::BeginRunSubsystem(const std::pair<SubsysReco *, PHCompositeNode *> &subsys)
{
  int iret = 0;
  std::ostringstream newdirname;
  newdirname << subsys.second->getName() << "/" << subsys.first->Name();
  if (!gROOT->cd(newdirname.str().c_str()))
  {
    std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
              << subsys.second->getName()
              << " - send e-mail to off-l with your macro" << std::endl;
    exit(1);
  }
  else
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      std::cout << "BeginRun: cded to " << newdirname.str().c_str() << std::endl;
    }
  }

  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Fun4AllServer::BeginRun: InitRun for " << subsys.first->Name() << std::endl;
  }
  try
  {
#ifdef FFAMEMTRACKER
    ffamemtracker->Start(subsys.first->Name(), "SubsysReco");
#endif
    iret = subsys.first->InitRun(subsys.second);
#ifdef FFAMEMTRACKER
    ffamemtracker->Stop(subsys.first->Name(), "SubsysReco");
#endif
  }
  catch (const std::exception &e)
  {
    std::cout << PHWHERE << " caught exception thrown during SubsysReco::InitRun() from "
              << subsys.first->Name() << std::endl;
    std::cout << "error: " << e.what() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::InitRun() from "
              << subsys.first->Name() << std::endl;
    exit(1);
  }

  if (iret == Fun4AllReturnCodes::ABORTRUN)
  {
    std::cout << PHWHERE << "Module " << subsys.first->Name() << " issued Abort Run, exiting" << std::endl;
    exit(-1);
  }
  else if (iret != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << PHWHERE << "Module " << subsys.first->Name() << " issued non Fun4AllReturnCodes::EVENT_OK return code " << iret << " in InitRun()" << std::endl;
    exit(-2);
  }
  return iret;
}

int Fun4AllServer::CountOutNodes(PHCompositeNode *startNode)
{
  int icount = 0;
  icount = CountOutNodesRecursive(startNode, icount);
  return icount;
}

int Fun4AllServer::CountOutNodesRecursive(PHCompositeNode *startNode, const int icount)
{
  PHNodeIterator nodeiter(startNode);
  PHPointerListIterator<PHNode> iterat(nodeiter.ls());
  PHNode *thisNode;
  int icnt = icount;
  while ((thisNode = iterat()))
  {
    if ((thisNode->getType() == "PHCompositeNode"))
    {
      icnt = CountOutNodesRecursive(static_cast<PHCompositeNode *>(thisNode), icnt);  // if this is a CompositeNode do this trick again
    }
    else
    {
      icnt++;
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << thisNode->getName() << ", Node Count: " << icnt << std::endl;
      }
    }
  }
  return icnt;
}

int Fun4AllServer::MakeNodesTransient(PHCompositeNode *startNode)
{
  PHNodeIterator nodeiter(startNode);
  PHPointerListIterator<PHNode> iterat(nodeiter.ls());
  PHNode *thisNode;
  while ((thisNode = iterat()))
  {
    if ((thisNode->getType() == "PHCompositeNode"))
    {
      MakeNodesTransient(static_cast<PHCompositeNode *>(thisNode));  // if this is a CompositeNode do this trick again
    }
    else
    {
      thisNode->makeTransient();
    }
  }
  return 0;
}

int Fun4AllServer::MakeNodesPersistent(PHCompositeNode *startNode)
{
  PHNodeIterator nodeiter(startNode);
  PHPointerListIterator<PHNode> iterat(nodeiter.ls());
  PHNode *thisNode;
  while ((thisNode = iterat()))
  {
    if ((thisNode->getType() == "PHCompositeNode"))
    {
      MakeNodesPersistent(static_cast<PHCompositeNode *>(thisNode));  // if this is a CompositeNode do this trick again
    }
    else
    {
      thisNode->makePersistent();
    }
  }
  return 0;
}

int Fun4AllServer::EndRun(const int runno)
{
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::iterator iter;
  gROOT->cd(default_Tdirectory.c_str());
  std::string currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Fun4AllServer::EndRun: EndRun("
                << runno << ") for " << (*iter).first->Name() << std::endl;
    }
    std::ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
                << (*iter).second->getName()
                << " - send e-mail to off-l with your macro" << std::endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "EndRun: cded to " << newdirname.str().c_str() << std::endl;
      }
    }
    try
    {
      (*iter).first->EndRun(runno);
    }
    catch (const std::exception &e)
    {
      std::cout << PHWHERE << " caught exception thrown during SubsysReco::EndRun() from "
                << (*iter).first->Name() << std::endl;
      std::cout << "error: " << e.what() << std::endl;
      exit(1);
    }
    catch (...)
    {
      std::cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::EndRun() from "
                << (*iter).first->Name() << std::endl;
      exit(1);
    }
  }
  gROOT->cd(currdir.c_str());

  return 0;
}

int Fun4AllServer::End()
{
  recoConsts *rc = recoConsts::instance();
  EndRun(rc->get_IntFlag("RUNNUMBER"));  // call SubsysReco EndRun methods for current run
  int i = 0;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::iterator iter;
  gROOT->cd(default_Tdirectory.c_str());
  std::string currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      std::cout << "Fun4AllServer::End: End for " << (*iter).first->Name() << std::endl;
    }
    std::ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      std::cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
                << (*iter).second->getName()
                << " - send e-mail to off-l with your macro" << std::endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        std::cout << "End: cded to " << newdirname.str().c_str() << std::endl;
      }
    }
    try
    {
      i += (*iter).first->End((*iter).second);
    }
    catch (const std::exception &e)
    {
      std::cout << PHWHERE << " caught exception thrown during SusbsysReco::End() from "
                << (*iter).first->Name() << std::endl;
      std::cout << "error: " << e.what() << std::endl;
      exit(1);
    }
    catch (...)
    {
      std::cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::End() from "
                << (*iter).first->Name() << std::endl;
      exit(1);
    }
  }
  gROOT->cd(currdir.c_str());
  PHNodeIterator nodeiter(TopNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(nodeiter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << "No Run Node, not writing Runwise info" << std::endl;
  }
  else
  {
    if (!OutputManager.empty())  // there are registered IO managers
    {
      MakeNodesTransient(runNode);  // make all nodes transient by default
      std::vector<Fun4AllOutputManager *>::iterator IOiter;
      for (IOiter = OutputManager.begin(); IOiter != OutputManager.end(); ++IOiter)
      {
        (*IOiter)->WriteNode(runNode);
      }
    }
  }
  // close output files (check for existing output managers is
  // done inside outfileclose())
  outfileclose();

  if (ScreamEveryEvent)
  {
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "Now that we are at the End(), please fix the following "
              << ScreamEveryEvent << " problem(s):" << std::endl;
    std::vector<std::string>::const_iterator viter;
    for (viter = ComplaintList.begin(); viter != ComplaintList.end(); ++viter)
    {
      std::cout << *viter << std::endl;
    }
    std::cout << " " << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
    std::cout << "*******************************************************************************" << std::endl;
  }

  return i;
}

void Fun4AllServer::Print(const std::string &what) const
{
  if (what == "ALL" || what == "HISTOS")
  {
    // loop over the map and print out the content (name and location in memory)
    for (auto &histoman : HistoManager)
    {
      histoman->Print(what);
    }
  }
  if (what == "ALL" || what == "SUBSYSTEMS")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Subsystems in Fun4AllServer:" << std::endl;

    std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::const_iterator miter;
    for (miter = Subsystems.begin(); miter != Subsystems.end(); ++miter)
    {
      std::cout << (*miter).first->Name()
                << " running under topNode " << (*miter).second->getName() << std::endl;
    }
    std::cout << std::endl;
  }

  if (what == "ALL" || what == "INPUTMANAGER")
  {
    // the input managers are managed by the input singleton
    for (auto &syncman : SyncManagers)
    {
      std::cout << "SyncManager: " << syncman->Name() << std::endl;
      syncman->Print(what);
    }
  }

  if (what == "ALL" || what.find("OUTPUTMANAGER") != std::string::npos)
  {
    // loop over the map and print out the content (name and location in memory)
    std::string pass_on = what;
    if (pass_on == "ALL" || pass_on == "OUTPUTMANAGER")
    {
      std::cout << "--------------------------------------" << std::endl
                << std::endl;
      std::cout << "List of OutputManagers in Fun4AllServer:" << std::endl;
      pass_on = "ALL";
    }
    else
    {
      std::string::size_type pos = pass_on.find('%');
      pass_on = pass_on.substr(pos + 1, pass_on.size());
    }
    for (auto &outman : OutputManager)
    {
      outman->Print(pass_on);
    }
    std::cout << std::endl;
  }
  if (what == "ALL" || what == "TOPNODES")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of TopNodes in Fun4AllServer:" << std::endl;

    std::map<std::string, PHCompositeNode *>::const_iterator iter;
    for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
    {
      std::cout << iter->first << " is at " << std::hex
                << iter->second << std::dec << std::endl;
    }
    std::cout << std::endl;
  }
  if (what == "ALL" || what == "NODETREE")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of Nodes in Fun4AllServer:" << std::endl;

    std::map<std::string, PHCompositeNode *>::const_iterator iter;
    for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
    {
      std::cout << "Node Tree under TopNode " << iter->first << std::endl;
      PHNodeIterator nodeiter(iter->second);
      nodeiter.print();
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  return;
}

void Fun4AllServer::identify(std::ostream &out) const
{
  out << "Fun4AllServer Name: " << Name() << std::endl;
  return;
}

int Fun4AllServer::outfileclose()
{
  while (!OutputManager.empty())
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Erasing OutputManager "
                << (*OutputManager.begin())->Name()
                << " at memory location " << *(OutputManager.begin()) << std::endl;
    }
    delete *(OutputManager.begin());
    OutputManager.erase(OutputManager.begin());
  }
  return 0;
}

int Fun4AllServer::InitNodeTree(PHCompositeNode *topNode)
{
  PHCompositeNode *dstNode, *runNode, *parNode;
  dstNode = new PHCompositeNode("DST");
  topNode->addNode(dstNode);
  runNode = new PHCompositeNode("RUN");
  topNode->addNode(runNode);
  parNode = new PHCompositeNode("PAR");
  topNode->addNode(parNode);
  return 0;
}

PHCompositeNode *
Fun4AllServer::topNode(const std::string &name)
{
  std::map<std::string, PHCompositeNode *>::const_iterator iter;
  iter = topnodemap.find(name);
  if (iter != topnodemap.end())
  {
    return iter->second;
  }
  AddTopNode(name);
  iter = topnodemap.find(name);
  if (iter != topnodemap.end())
  {
    InitNodeTree(iter->second);
    return iter->second;
  }
  std::cout << PHWHERE << " Could not create new topNode " << name
            << " send email to off-l with the following printout: " << std::endl;
  for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
  {
    std::cout << iter->first << " is at " << std::hex << iter->second << std::dec << std::endl;
  }
  exit(1);
}

int Fun4AllServer::AddTopNode(const std::string &name)
{
  std::map<std::string, PHCompositeNode *>::const_iterator iter;
  iter = topnodemap.find(name);
  if (iter != topnodemap.end())
  {
    return -1;
  }
  PHCompositeNode *newNode = new PHCompositeNode(name.c_str());
  topnodemap[name] = newNode;
  return 0;
}

PHCompositeNode *Fun4AllServer::getNode(const std::string &name, const std::string &topnodename)
{
  PHNodeIterator iter(topNode(topnodename));
  PHCompositeNode *thisNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", name));
  if (!thisNode)
  {
    thisNode = new PHCompositeNode(name);
    topNode(topnodename)->addNode(thisNode);
  }
  return thisNode;
}

int Fun4AllServer::registerInputManager(Fun4AllInputManager *InManager)
{
  int iret = defaultSyncManager->registerInputManager(InManager);
  return iret;
}

Fun4AllInputManager *
Fun4AllServer::getInputManager(const std::string &name)
{
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    if ((*iter)->getInputManager(name))
    {
      return (*iter)->getInputManager(name);
    }
  }
  std::cout << "Could not locate input manager " << name << std::endl;
  return nullptr;
}

int Fun4AllServer::PrdfEvents() const
{
  return (defaultSyncManager->PrdfEvents());
}

int Fun4AllServer::DstEvents() const
{
  return (defaultSyncManager->TotalEvents());
}

//_________________________________________________________________
int Fun4AllServer::run(const int nevnts, const bool require_nevents)
{
  recoConsts *rc = recoConsts::instance();
  static bool run_number_forced = rc->FlagExist("RUNNUMBER");
  static int ifirst = 1;
  if (ifirst && run_number_forced)
  {
    runnumber = rc->get_IntFlag("RUNNUMBER");
    std::cout << "Fun4AllServer: Runnumber forced to " << runnumber << " by RUNNUMBER IntFlag" << std::endl;
  }
  int iret = 0;
  int icnt = 0;
  int icnt_good = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  while (!iret)
  {
    int resetnodetree = 0;
    for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
    {
      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << "executing run for input master " << (*iter)->Name() << std::endl;
      }
      int retval = (*iter)->run(1);
      // if a new input file is opened during syncing and it contains
      // different nodes
      // as the previous one, the info in the nodes which are only in
      // the previous file will be carried to another event. We also
      // do not know under which topNode the input managers put
      // their data. This is why
      // the whole node tree is resetted whenever one of the Sync Managers
      // requires it.
      if (retval == Fun4AllReturnCodes::RESET_NODE_TREE)
      {
        resetnodetree = 1;
      }
      else
      {
        iret += retval;
      }
    }
    if (resetnodetree)
    {
      // if the node tree needs resetting, we just push the current
      // event(s) (which are all properly synced at this point)
      // back into the input managers and just read again.
      for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
      {
        (*iter)->PushBackInputMgrsEvents(1);
      }
      ResetNodeTree();
      continue;  // go back to run loop
    }
    if (iret)
    {
      break;
    }
    int currentrun = 0;
    for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
    {
      int runno = (*iter)->CurrentRun();
      //	  std::cout << (*iter)->Name() << " run no: " << runno << std::endl;
      if (runno != 0)
      {
        if (currentrun == 0)
        {
          currentrun = runno;
        }
        else
        {
          if (currentrun != runno)
          {
            std::cout << "Mixing of Runs within same event is not supported" << std::endl;
            std::cout << "Here is the list of Sync Managers and their runnumbers:" << std::endl;
            std::vector<Fun4AllSyncManager *>::const_iterator syiter;
            for (syiter = SyncManagers.begin(); syiter != SyncManagers.end(); ++syiter)
            {
              std::cout << (*syiter)->Name() << " run number: " << (*syiter)->CurrentRun() << std::endl;
            }
            std::cout << "Exiting now" << std::endl;
            exit(1);
          }
        }
      }
    }
    if (ifirst)
    {
      if (currentrun != runnumber && !run_number_forced)  // use real run if not forced
      {
        runnumber = currentrun;
      }
      setRun(runnumber);
      BeginRun(runnumber);
      ifirst = 0;
    }
    else if (!run_number_forced)
    {
      if (currentrun != runnumber)
      {
        EndRun(runnumber);
        runnumber = currentrun;
        setRun(runnumber);
        BeginRun(runnumber);
      }
    }
    if (Verbosity() >= 1)
    {
      std::cout << "Fun4AllServer::run - processing event "
                << (icnt + 1) << " from run " << runnumber << std::endl;
    }

    if (icnt == 0 and Verbosity() > VERBOSITY_QUIET)
    {
      // increase verbosity for the first event in verbose modes
      int iverb = Verbosity();
      Verbosity(++iverb);
    }

    iret = process_event();

    if (icnt == 0 and Verbosity() > VERBOSITY_QUIET)
    {
      // increase verbosity for the first event in verbose modes
      int iverb = Verbosity();
      Verbosity(--iverb);
    }

    ++icnt;  // completed one event processing

    if (require_nevents)
    {
      if (std::find(RetCodes.begin(),
                    RetCodes.end(),
                    static_cast<int>(Fun4AllReturnCodes::ABORTEVENT)) == RetCodes.end())
        icnt_good++;
      if (iret || (nevnts > 0 && icnt_good >= nevnts))
        break;
    }
    else if (iret || (nevnts > 0 && icnt >= nevnts))
    {
      break;
    }
  }
  return iret;
}

//_________________________________________________________________
int Fun4AllServer::skip(const int nevnts)
{
  int iret = 0;
  if (nevnts > 0)  // do not execute for nevnts <= 0
  {
    std::vector<Fun4AllSyncManager *>::const_iterator iter;
    for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
    {
      iret += (*iter)->skip(nevnts);
    }
    eventcounter += nevnts;  // update event counter so it reflects the number of events in the input
  }
  return iret;
}

//_________________________________________________________________
int Fun4AllServer::fileopen(const std::string &managername, const std::string &filename)
{
  int iret = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->fileopen(managername, filename);
  }
  return iret;
}

int Fun4AllServer::BranchSelect(const std::string &managername, const std::string &branch, int iflag)
{
  int iret = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->BranchSelect(managername, branch, iflag);
  }
  return iret;
}

int Fun4AllServer::BranchSelect(const std::string &branch, int iflag)
{
  int iret = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->BranchSelect(branch, iflag);
  }
  return iret;
}

int Fun4AllServer::setBranches(const std::string &managername)
{
  int iret = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->setBranches(managername);
  }
  return iret;
}

int Fun4AllServer::setBranches()
{
  int iret = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->setBranches();
  }
  return iret;
}

int Fun4AllServer::fileclose(const std::string &managername)
{
  int iret = 0;
  std::vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->fileclose(managername);
  }
  return iret;
}

int Fun4AllServer::SegmentNumber()
{
  int iret = defaultSyncManager->SegmentNumber();
  return iret;
}

void Fun4AllServer::GetInputFullFileList(std::vector<std::string> &fnames) const
{
  defaultSyncManager->GetInputFullFileList(fnames);
  return;
}

int Fun4AllServer::DisconnectDB()
{
  return 0;
}

unsigned
Fun4AllServer::GetTopNodes(std::vector<std::string> &names) const
{
  std::map<std::string, PHCompositeNode *>::const_iterator iter;
  for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
  {
    names.push_back(iter->first);
  }
  return names.size();
}

void Fun4AllServer::GetOutputManagerList(std::vector<std::string> &names) const
{
  names.clear();
  std::vector<Fun4AllOutputManager *>::const_iterator iter;
  for (iter = OutputManager.begin(); iter != OutputManager.end(); ++iter)
  {
    names.push_back((*iter)->Name());
  }
  return;
}

void Fun4AllServer::GetModuleList(std::vector<std::string> &names) const
{
  names.clear();
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>>::const_iterator iter;
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    names.push_back((*iter).first->Name());
  }
  return;
}

int Fun4AllServer::registerSyncManager(Fun4AllSyncManager *newmaster)
{
  for (auto &syncman : SyncManagers)
  {
    if (syncman->Name() == newmaster->Name())
    {
      std::cout << "Input Master " << newmaster->Name()
                << " already registered" << std::endl;
      return -1;
    }
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Registering Input Master " << newmaster->Name() << std::endl;
  }
  SyncManagers.push_back(newmaster);
  return 0;
}

Fun4AllSyncManager *
Fun4AllServer::getSyncManager(const std::string &name)
{
  if (name == "DefaultSyncManager")
  {
    return defaultSyncManager;
  }
  std::vector<Fun4AllSyncManager *>::iterator iter;

  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    if ((*iter)->Name() == name)
    {
      return *iter;
    }
  }
  std::cout << "Could not find Input Master " << name << std::endl;
  return nullptr;
}

int Fun4AllServer::setRun(const int runno)
{
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runno);
  PHTimeStamp *tstamp = nullptr;
  if (!tstamp)
  {
    tstamp = new PHTimeStamp(0);
    std::cout << "Fun4AllServer::setRun(): could not get timestamp for run  " << runno
              << ", using tics(0) timestamp: ";
    tstamp->print();
    std::cout << std::endl;
  }
  delete tstamp;
  FrameWorkVars->SetBinContent(RUNNUMBERBIN, (Stat_t) runno);
  return 0;
}

void Fun4AllServer::NodeIdentify(const std::string &name)
{
  PHObject *obj = findNode::getClass<PHObject>(TopNode, name);
  if (obj)
  {
    obj->identify();
  }
  else
  {
    std::cout << "Could not locate node " << name
              << " or no PHObject Node" << std::endl;
  }
  return;
}

void Fun4AllServer::PrintTimer(const std::string &name)
{
  std::map<const std::string, PHTimer>::const_iterator iter;
  PHTimer::PRINT(std::cout, "**");
  if (name.empty())
  {
    for (iter = timer_map.begin(); iter != timer_map.end(); ++iter)
    {
      iter->second.print_stat();
    }
  }
  else
  {
    iter = timer_map.find(name);
    if (iter != timer_map.end())
    {
      iter->second.print_stat();
    }
    else
    {
      std::cout << "No timer with name " << name << " found" << std::endl;
      std::cout << "Existing timers:" << std::endl;
      for (iter = timer_map.begin(); iter != timer_map.end(); ++iter)
      {
        std::cout << iter->first << std::endl;
      }
    }
  }
  return;
}

void Fun4AllServer::PrintMemoryTracker(const std::string &name) const
{
#ifdef FFAMEMTRACKER
  ffamemtracker->PrintMemoryTracker(name);
#else
  std::cout << "PrintMemoryTracker called with " << name << " is disabled" << std::endl;
#endif
  return;
}
