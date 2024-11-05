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
      std::cout << PHWHERE << " Error creating TDirectory topdir " << topnodename << std::endl;
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
      std::cout << PHWHERE << "Error creating TDirectory subdir " << subsystem->Name() << std::endl;
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
    if (Verbosity() >= 3)
    {
      std::cout << "Calling Init() for Subsystem " << subsystem->Name() << std::endl;
    }
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
      // NOLINTNEXTLINE(hicpp-avoid-goto)
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
  for (auto &Subsystem : Subsystems)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Fun4AllServer::process_event processing " << Subsystem.first->Name() << std::endl;
    }
    std::string newdirname = Subsystem.second->getName() + "/" + Subsystem.first->Name();
    if (!gROOT->cd(newdirname.c_str()))
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
        std::cout << "process_event: cded to " << newdirname << std::endl;
      }
    }

    PHTimer subsystem_timer("SubsystemTimer");
    subsystem_timer.restart();

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
      else if (RetCodes[icnt] == Fun4AllReturnCodes::ABORTPROCESSING)
      {
        eventbad = 1;
        retcodesmap[Fun4AllReturnCodes::ABORTPROCESSING]++;
        std::cout << "Fun4AllServer::Abort Processing by " << Subsystem.first->Name() << std::endl;
        return Fun4AllReturnCodes::ABORTPROCESSING;
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
    subsystem_timer.stop();
    double TimeSubsystem = subsystem_timer.elapsed();
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Fun4AllServer::process_event processing " << Subsystem.first->Name()
                << " processing total time: " << TimeSubsystem << " ms" << std::endl;
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
          if ((*iterOutMan)->EventsWritten() >= (*iterOutMan)->GetNEvents())
          {
            if (Verbosity() > 0)
            {
              std::cout << PHWHERE << (*iterOutMan)->Name() << " wrote " << (*iterOutMan)->EventsWritten()
                        << " events, closing " << (*iterOutMan)->OutFileName() << std::endl;
            }
            PHNodeIterator nodeiter(TopNode);
            PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(nodeiter.findFirst("PHCompositeNode", "RUN"));
            MakeNodesTransient(runNode);  // make all nodes transient by default
            (*iterOutMan)->WriteNode(runNode);
            (*iterOutMan)->RunAfterClosing();
          }
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
  for (auto &Subsystem : Subsystems)
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
    for (const auto &nodename : ResetNodeList)
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
    // rc->set_TimeStamp(*beginruntimestamp);
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
  std::string newdirname = subsys.second->getName() + "/" + subsys.first->Name();
  if (!gROOT->cd(newdirname.c_str()))
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
      std::cout << "BeginRun: cded to " << newdirname << std::endl;
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

// NOLINTNEXTLINE(misc-no-recursion)
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

// NOLINTNEXTLINE(misc-no-recursion)
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

// NOLINTNEXTLINE(misc-no-recursion)
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
    std::string newdirname = (*iter).second->getName() + "/" + (*iter).first->Name();
    if (!gROOT->cd(newdirname.c_str()))
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
        std::cout << "EndRun: cded to " << newdirname << std::endl;
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
  recoConsts *rc = recoCo