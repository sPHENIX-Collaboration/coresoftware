#include "Fun4AllServer.h"
#include "Fun4AllHistoBinDefs.h"
#include "Fun4AllInputManager.h"
#include "Fun4AllOutputManager.h"
#include "Fun4AllReturnCodes.h"
#include "Fun4AllSyncManager.h"
#include "SubsysReco.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIOManager.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHNodeReset.h>
#include <phool/PHObject.h>
#include <phool/PHPointerListIterator.h>
#include <phool/PHTimeStamp.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include <boost/foreach.hpp>

#include <sys/utsname.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>

using namespace std;

Fun4AllServer *Fun4AllServer::__instance = 0;

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
  , OutNodeCount(0)
  , bortime_override(0)
  , ScreamEveryEvent(0)
  , unregistersubsystem(0)
  , runnumber(0)
  , eventnumber(0)
  , beginruntimestamp(nullptr)
  , keep_db_connected(0)
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
  ostringstream histomanagername;
  histomanagername << Name() << "HISTOS";
  ServerHistoManager = new Fun4AllHistoManager(histomanagername.str());
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

int Fun4AllServer::dumpHistos(const string &filename, const string &openmode)
{
  int iret = 0;
  cout << "Fun4AllServer::dumpHistos() dumping histograms" << endl;
  if (!filename.empty())
  {
    ServerHistoManager->setOutfileName(filename);
  }
  vector<Fun4AllHistoManager *>::const_iterator hiter;
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

bool Fun4AllServer::registerHisto(const char *hname, TNamed *h1d, const int replace)
{
  return ServerHistoManager->registerHisto(hname, h1d, replace);
}

int Fun4AllServer::isHistoRegistered(const string &name) const
{
  int iret = ServerHistoManager->isHistoRegistered(name);
  return iret;
}

int Fun4AllServer::registerSubsystem(SubsysReco *subsystem, const string &topnodename)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  // if somebody opens a TFile (or changes the gDirectory) in the ctor
  // we need to set it to a "known" directory
  gROOT->cd(default_Tdirectory.c_str());
  string currdir = gDirectory->GetPath();
  TDirectory *tmpdir = gDirectory;
  if (!tmpdir->FindObject(topnodename.c_str()))
  {
    tmpdir = tmpdir->mkdir(topnodename.c_str());
    if (!tmpdir)
    {
      cout << "Error creating TDirectory topdir " << topnodename.c_str() << endl;
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
      cout << "Error creating TDirectory subdir " << subsystem->Name() << endl;
      exit(1);
    }
    // store the TDir pointer so it can be cleaned up in the dtor
    // if one deletes it here the Histograms are dangling somewhere
    // in root
    TDirCollection.push_back(tmpdir);
  }
  PHCompositeNode *subsystopNode = se->topNode(topnodename);
  pair<SubsysReco *, PHCompositeNode *> newsubsyspair(subsystem, subsystopNode);
  int iret = 0;
  try
  {
    iret = subsystem->Init(subsystopNode);
  }
  catch (const exception &e)
  {
    cout << PHWHERE << " caught exception thrown during SubsysReco::Init() from "
         << subsystem->Name() << endl;
    cout << "error: " << e.what() << endl;
    exit(1);
  }
  catch (...)
  {
    cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::Init() from "
         << subsystem->Name() << endl;
    exit(1);
  }
  gROOT->cd(currdir.c_str());
  if (iret)
  {
    if (iret == Fun4AllReturnCodes::DONOTREGISTERSUBSYSTEM)
    {
      if (Verbosity() >= VERBOSITY_SOME)
      {
        cout << "Not Registering Subsystem " << subsystem->Name() << endl;
      }
      return 0;
    }
    cout << PHWHERE << " Error initializing subsystem "
         << subsystem->Name() << ", return code: " << iret << endl;
    return iret;
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "Registering Subsystem " << subsystem->Name() << endl;
  }
  Subsystems.push_back(newsubsyspair);
  ostringstream timer_name;
  timer_name << subsystem->Name() << "_" << topnodename;
  PHTimer timer(timer_name.str());
  if (timer_map.find(timer_name.str()) == timer_map.end())
  {
    timer_map[timer_name.str()] = timer;
  }
  RetCodes.push_back(iret);  // vector with return codes
  return 0;
}

int Fun4AllServer::unregisterSubsystem(SubsysReco *subsystem)
{
  pair<SubsysReco *, PHCompositeNode *> subsyspair(subsystem, 0);
  DeleteSubsystems.push_back(subsyspair);
  unregistersubsystem = 1;
  return 0;
}

int Fun4AllServer::unregisterSubsystemsNow()
{
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator sysiter, removeiter;
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
      cout << "unregisterSubsystem: Could not find SubsysReco "
           << (*removeiter).first->Name()
           << " in Fun4All Reco Module list" << endl;
      delete (*removeiter).first;
      continue;
    }
    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "Removing Subsystem: " << (*removeiter).first->Name()
           << " at index " << index << endl;
    }
    Subsystems.erase(Subsystems.begin() + index);
    delete (*removeiter).first;
    // also update the vector with return codes
    RetCodes.erase(RetCodes.begin() + index);
    vector<Fun4AllOutputManager *>::iterator outiter;
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
Fun4AllServer::getSubsysReco(const string &name)
{
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator sysiter;
  for (sysiter = Subsystems.begin(); sysiter != Subsystems.end(); ++sysiter)
  {
    if ((*sysiter).first->Name() == name)
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "Found Subsystem " << name << endl;
      }
      return (*sysiter).first;
    }
  }
  cout << "Could not find SubsysReco " << name << endl;
  return 0;
}

int Fun4AllServer::AddComplaint(const string &complaint, const string &remedy)
{
  ScreamEveryEvent++;
  string separatorstring = "------------------------------";
  ostringstream complaintno;
  complaintno << "Problem No " << ScreamEveryEvent;

  ComplaintList.push_back(separatorstring);
  ComplaintList.push_back(complaintno.str());
  ComplaintList.push_back(complaint);
  ComplaintList.push_back(" ");
  ComplaintList.push_back("Remedy:");
  ComplaintList.push_back(remedy);
  ComplaintList.push_back(separatorstring);
  return 0;
}

int Fun4AllServer::registerOutputManager(Fun4AllOutputManager *manager)
{
  vector<Fun4AllOutputManager *>::iterator iter;
  for (iter = OutputManager.begin(); iter != OutputManager.end(); ++iter)
  {
    if ((*iter)->Name() == manager->Name())
    {
      cout << "OutputManager " << manager->Name() << " allready in list" << endl;
      return -1;
    }
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "Registering OutputManager " << manager->Name() << endl;
  }
  UpdateEventSelector(manager);
  OutputManager.push_back(manager);
  return 0;
}

int Fun4AllServer::UpdateEventSelector(Fun4AllOutputManager *manager)
{
  vector<string>::iterator striter;
  vector<pair<SubsysReco *, PHCompositeNode *> >::const_iterator subsysiter;

tryagain:
  manager->RecoModuleIndex()->clear();
  for (striter = manager->EventSelector()->begin(); striter != manager->EventSelector()->end(); ++striter)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      cout << PHWHERE << "striter: " << *striter << endl;
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
          cout << PHWHERE << "setting RecoModuleIndex to " << index << endl;
        }
        found = 1;
        break;
      }
      index++;
    }
    if (!found)
    {
      cout << "Could not find module " << *striter
           << ", removing it from list of event selector modules" << endl;
      manager->EventSelector()->erase(striter);
      goto tryagain;
    }
  }
  return 0;
}

Fun4AllOutputManager *
Fun4AllServer::getOutputManager(const string &name)
{
  vector<Fun4AllOutputManager *>::iterator iter;
  for (iter = OutputManager.begin(); iter != OutputManager.end(); ++iter)
  {
    if (name == (*iter)->Name())
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "Found OutputManager " << name << endl;
      }
      return *iter;
    }
  }
  cout << "Could not find OutputManager" << name << endl;
  return 0;
}

Fun4AllHistoManager *
Fun4AllServer::getHistoManager(const string &name)
{
  vector<Fun4AllHistoManager *>::iterator iter;
  for (iter = HistoManager.begin(); iter != HistoManager.end(); ++iter)
  {
    if ((*iter)->Name() == name)
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "Found HistoManager " << name << endl;
      }
      return *iter;
    }
  }
  if (Verbosity() >= VERBOSITY_MORE)
  {
    cout << "Could not find HistoManager " << name << endl;
  }
  return 0;
}

int Fun4AllServer::registerHistoManager(Fun4AllHistoManager *manager)
{
  vector<Fun4AllHistoManager *>::iterator iter;
  for (iter = HistoManager.begin(); iter != HistoManager.end(); ++iter)
  {
    if ((*iter)->Name() == manager->Name())
    {
      cout << "HistoManager " << manager->Name() << " allready in list" << endl;
      return -1;
    }
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "Registering HistoManager " << manager->Name() << endl;
  }
  HistoManager.push_back(manager);
  return 0;
}

TNamed *
Fun4AllServer::getHisto(const unsigned int ihisto) const
{
  return ServerHistoManager->getHisto(ihisto);
}

const char *
Fun4AllServer::getHistoName(const unsigned int ihisto) const
{
  return (ServerHistoManager->getHistoName(ihisto));
}

TNamed *Fun4AllServer::getHisto(const string &hname) const
{
  return (ServerHistoManager->getHisto(hname));
}

int Fun4AllServer::process_event(PHCompositeNode * /*topNode*/)
{
  int iret = process_event();
  return iret;
}

int Fun4AllServer::process_event()
{
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  unsigned icnt = 0;
  int eventbad = 0;
  if (ScreamEveryEvent)
  {
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
    cout << "Now that I have your attention, please fix the following "
         << ScreamEveryEvent << " problem(s):" << endl;
    vector<string>::const_iterator viter;
    for (viter = ComplaintList.begin(); viter != ComplaintList.end(); ++viter)
    {
      cout << *viter << endl;
    }
    cout << " " << endl;
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
  }
  if (unregistersubsystem)
  {
    unregisterSubsystemsNow();
  }
  gROOT->cd(default_Tdirectory.c_str());
  string currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "Fun4AllServer::process_event processing " << (*iter).first->Name() << endl;
    }
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "process_event: cded to " << newdirname.str().c_str() << endl;
      }
    }

    try
    {
      ostringstream timer_name;
      timer_name << (*iter).first->Name() << "_" << (*iter).second->getName();
      std::map<const std::string, PHTimer>::iterator titer = timer_map.find(timer_name.str());
      bool timer_found = false;
      if (titer != timer_map.end())
      {
        timer_found = true;
	titer->second.restart();
      }
      else
      {
	cout << "could not find timer for " << timer_name << endl;
      }
      RetCodes[icnt] = (*iter).first->process_event((*iter).second);
      if (timer_found)
      {
	titer->second.stop();
      }

    }
    catch (const exception &e)
    {
      cout << PHWHERE << " caught exception thrown during process_event from "
           << (*iter).first->Name() << endl;
      cout << "error: " << e.what() << endl;
      exit(1);
    }
    catch (...)
    {
      cout << PHWHERE << " caught unknown type exception thrown during process_event from "
           << (*iter).first->Name() << endl;
      exit(1);
    }
    if (RetCodes[icnt])
    {
      if (RetCodes[icnt] == Fun4AllReturnCodes::DISCARDEVENT)
      {
        if (Verbosity() >= VERBOSITY_EVEN_MORE)
        {
          cout << "Fun4AllServer::Discard Event by " << (*iter).first->Name() << endl;
        }
      }
      else if (RetCodes[icnt] == Fun4AllReturnCodes::ABORTEVENT)
      {
        retcodesmap[Fun4AllReturnCodes::ABORTEVENT]++;
        eventbad = 1;
        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "Fun4AllServer::Abort Event by " << (*iter).first->Name() << endl;
        }
        break;
      }
      else if (RetCodes[icnt] == Fun4AllReturnCodes::ABORTRUN)
      {
        retcodesmap[Fun4AllReturnCodes::ABORTRUN]++;
        cout << "Fun4AllServer::Abort Run by " << (*iter).first->Name() << endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      else
      {
        cout << "Fun4AllServer::Unknown return code: "
             << RetCodes[icnt] << " from process_event method of "
             << (*iter).first->Name() << endl;
        cout << "This smells like an uninitialized return code and" << endl;
        cout << "it is too dangerous to continue, this Run will be aborted" << endl;
        cout << "If you do not know how to fix this please send mail to" << endl;
        cout << "phenix-off-l with this message" << endl;
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
        cout << PHWHERE << " FATAL: Someone changed the number of Output Nodes on the fly, from " << OutNodeCount << " to " << newcount << endl;
        exit(1);
      }
      vector<Fun4AllOutputManager *>::iterator iterOutMan;
      for (iterOutMan = OutputManager.begin(); iterOutMan != OutputManager.end(); ++iterOutMan)
      {
        if (!(*iterOutMan)->DoNotWriteEvent(&RetCodes))
        {
          if (Verbosity() >= VERBOSITY_MORE)
          {
            cout << "Writing Event for " << (*iterOutMan)->Name() << endl;
          }
          (*iterOutMan)->WriteGeneric(dstNode);
        }
        else
        {
          if (Verbosity() >= VERBOSITY_MORE)
          {
            cout << "Not Writing Event for " << (*iterOutMan)->Name() << endl;
          }
        }
      }
    }
  }
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      cout << "Fun4AllServer::process_event Resetting Event " << (*iter).first->Name() << endl;
    }
    (*iter).first->ResetEvent((*iter).second);
  }
  BOOST_FOREACH (Fun4AllSyncManager *syncman, SyncManagers)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      cout << "Fun4AllServer::process_event Resetting Event for Sync Manager " << syncman->Name() << endl;
    }
    syncman->ResetEvent();
  }
  ResetNodeTree();
  return 0;
}

int Fun4AllServer::ResetNodeTree()
{
  vector<string> ResetNodeList;
  ResetNodeList.push_back("DST");
  PHNodeReset reset;
  reset.Verbosity(Verbosity()>0 ? Verbosity()-1: 0); // one lower verbosity level than Fun4AllServer
  map<string, PHCompositeNode *>::const_iterator iter;
  for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
  {
    PHNodeIterator mainIter((*iter).second);
    for (vector<string>::const_iterator nodename = ResetNodeList.begin();
         nodename != ResetNodeList.end(); ++nodename)
    {
      if (mainIter.cd(*nodename))
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
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_EVEN_MORE)
    {
      cout << "Fun4AllServer::Reset Resetting " << (*iter).first->Name() << endl;
    }
    i += (*iter).first->Reset((*iter).second);
  }
  vector<Fun4AllHistoManager *>::iterator hiter;
  for (hiter = HistoManager.begin(); hiter != HistoManager.end(); ++hiter)
  {
    (*hiter)->Reset();
  }
  return i;
}

int Fun4AllServer::BeginRunTimeStamp(PHTimeStamp &TimeStp)
{
  beginruntimestamp = new PHTimeStamp(TimeStp);
  cout << "Setting BOR timestamp to ";
  beginruntimestamp->print();
  cout << endl;
  bortime_override = 1;
  return 0;
}

int Fun4AllServer::BeginRun(const int runno)
{
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
    cout << "overriding BOR timestamp by ";
    beginruntimestamp->print();
    cout << endl;
    //rc->set_TimeStamp(*beginruntimestamp);
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "Fun4AllServer::BeginRun: Run number " << runno << " uses RECO TIMESTAMP: ";
    beginruntimestamp->print();
    cout << endl;
  }
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  int iret;

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
  string currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "BeginRun: cded to " << newdirname.str().c_str() << endl;
      }
    }

    if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "Fun4AllServer::BeginRun: InitRun for " << (*iter).first->Name() << endl;
    }
    try
    {
      iret = (*iter).first->InitRun((*iter).second);
    }
    catch (const exception &e)
    {
      cout << PHWHERE << " caught exception thrown during SubsysReco::InitRun() from "
           << (*iter).first->Name() << endl;
      cout << "error: " << e.what() << endl;
      exit(1);
    }
    catch (...)
    {
      cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::InitRun() from "
           << (*iter).first->Name() << endl;
      exit(1);
    }

    if (iret == Fun4AllReturnCodes::ABORTRUN)
    {
      cout << PHWHERE << "Module " << (*iter).first->Name() << " issued Abort Run, exiting" << endl;
      exit(-1);
    }
    else if (iret != Fun4AllReturnCodes::EVENT_OK)
    {
      cout << PHWHERE << "Module " << (*iter).first->Name() << " issued non Fun4AllReturnCodes::EVENT_OK return code " << iret << " in InitRun()" << endl;
      exit(-2);
    }
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
    cout << "WARNING WARNING, DBs will not be disconnected" << endl;
    cout << "This is for DB server testing purposes only" << endl;
    cout << "If you do not test our DB servers, remove" << endl;
    cout << "Fun4AllServer->KeepDBConnection()" << endl;
    cout << "from your macro" << endl;
  }
  // print out all node trees
  Print("NODETREE");
  return 0;
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
        cout << thisNode->getName() << ", Node Count: " << icnt << endl;
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
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  gROOT->cd(default_Tdirectory.c_str());
  string currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "Fun4AllServer::EndRun: EndRun("
           << runno << ") for " << (*iter).first->Name() << endl;
    }
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "EndRun: cded to " << newdirname.str().c_str() << endl;
      }
    }
    try
    {
      (*iter).first->EndRun(runno);
    }
    catch (const exception &e)
    {
      cout << PHWHERE << " caught exception thrown during SubsysReco::EndRun() from "
           << (*iter).first->Name() << endl;
      cout << "error: " << e.what() << endl;
      exit(1);
    }
    catch (...)
    {
      cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::EndRun() from "
           << (*iter).first->Name() << endl;
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
  vector<pair<SubsysReco *, PHCompositeNode *> >::iterator iter;
  gROOT->cd(default_Tdirectory.c_str());
  string currdir = gDirectory->GetPath();
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "Fun4AllServer::End: End for " << (*iter).first->Name() << endl;
    }
    ostringstream newdirname;
    newdirname << (*iter).second->getName() << "/" << (*iter).first->Name();
    if (!gROOT->cd(newdirname.str().c_str()))
    {
      cout << PHWHERE << "Unexpected TDirectory Problem cd'ing to "
           << (*iter).second->getName()
           << " - send e-mail to off-l with your macro" << endl;
      exit(1);
    }
    else
    {
      if (Verbosity() >= VERBOSITY_EVEN_MORE)
      {
        cout << "End: cded to " << newdirname.str().c_str() << endl;
      }
    }
    try
    {
      i += (*iter).first->End((*iter).second);
    }
    catch (const exception &e)
    {
      cout << PHWHERE << " caught exception thrown during SusbsysReco::End() from "
           << (*iter).first->Name() << endl;
      cout << "error: " << e.what() << endl;
      exit(1);
    }
    catch (...)
    {
      cout << PHWHERE << " caught unknown type exception thrown during SubsysReco::End() from "
           << (*iter).first->Name() << endl;
      exit(1);
    }
  }
  gROOT->cd(currdir.c_str());
  PHNodeIterator nodeiter(TopNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(nodeiter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    cout << "No Run Node, not writing Runwise info" << endl;
  }
  else
  {
    if (!OutputManager.empty())  // there are registered IO managers
    {
      vector<Fun4AllOutputManager *>::iterator IOiter;
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
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
    cout << "Now that we are at the End(), please fix the following "
         << ScreamEveryEvent << " problem(s):" << endl;
    vector<string>::const_iterator viter;
    for (viter = ComplaintList.begin(); viter != ComplaintList.end(); ++viter)
    {
      cout << *viter << endl;
    }
    cout << " " << endl;
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
    cout << "*******************************************************************************" << endl;
  }

  return i;
}

void Fun4AllServer::Print(const string &what) const
{
  if (what == "ALL" || what == "HISTOS")
  {
    // loop over the map and print out the content (name and location in memory)
    BOOST_FOREACH (Fun4AllHistoManager *histoman, HistoManager)
    {
      histoman->Print(what);
    }
  }
  if (what == "ALL" || what == "SUBSYSTEMS")
  {
    // loop over the map and print out the content (name and location in memory)
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of Subsystems in Fun4AllServer:" << endl;

    vector<pair<SubsysReco *, PHCompositeNode *> >::const_iterator miter;
    for (miter = Subsystems.begin(); miter != Subsystems.end(); ++miter)
    {
      cout << (*miter).first->Name()
           << " running under topNode " << (*miter).second->getName() << endl;
    }
    cout << endl;
  }

  if (what == "ALL" || what == "INPUTMANAGER")
  {
    // the input managers are managed by the input singleton
    BOOST_FOREACH (Fun4AllSyncManager *syncman, SyncManagers)
    {
      cout << "SyncManager: " << syncman->Name() << endl;
      syncman->Print(what);
    }
  }

  if (what == "ALL" || what.find("OUTPUTMANAGER") != string::npos)
  {
    // loop over the map and print out the content (name and location in memory)
    string pass_on = what;
    if (pass_on == "ALL" || pass_on == "OUTPUTMANAGER")
    {
      cout << "--------------------------------------" << endl
           << endl;
      cout << "List of OutputManagers in Fun4AllServer:" << endl;
      pass_on = "ALL";
    }
    else
    {
      string::size_type pos = pass_on.find("%");
      pass_on = pass_on.substr(pos + 1, pass_on.size());
    }
    BOOST_FOREACH (Fun4AllOutputManager *outman, OutputManager)
    {
      outman->Print(pass_on);
    }
    cout << endl;
  }
  if (what == "ALL" || what == "TOPNODES")
  {
    // loop over the map and print out the content (name and location in memory)
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of TopNodes in Fun4AllServer:" << endl;

    map<std::string, PHCompositeNode *>::const_iterator iter;
    for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
    {
      cout << iter->first << " is at " << hex
           << iter->second << dec << endl;
    }
    cout << endl;
  }
  if (what == "ALL" || what == "NODETREE")
  {
    // loop over the map and print out the content (name and location in memory)
    cout << "--------------------------------------" << endl
         << endl;
    cout << "List of Nodes in Fun4AllServer:" << endl;

    map<std::string, PHCompositeNode *>::const_iterator iter;
    for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
    {
      cout << "Node Tree under TopNode " << iter->first << endl;
      PHNodeIterator nodeiter(iter->second);
      nodeiter.print();
      cout << endl;
    }
    cout << endl;
  }
  return;
}

void Fun4AllServer::identify(std::ostream &out) const
{
  out << "Fun4AllServer Name: " << Name() << endl;
  return;
}

int Fun4AllServer::outfileclose()
{
  while (!OutputManager.empty())
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "Erasing OutputManager "
           << (*OutputManager.begin())->Name()
           << " at memory location " << *(OutputManager.begin()) << endl;
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
Fun4AllServer::topNode(const string &name)
{
  map<string, PHCompositeNode *>::const_iterator iter;
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
  cout << PHWHERE << " Could not create new topNode " << name
       << " send email to off-l with the following printout: " << endl;
  for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
  {
    cout << iter->first << " is at " << hex << iter->second << dec << endl;
  }
  exit(1);
}

int Fun4AllServer::AddTopNode(const string &name)
{
  map<string, PHCompositeNode *>::const_iterator iter;
  iter = topnodemap.find(name);
  if (iter != topnodemap.end())
  {
    return -1;
  }
  PHCompositeNode *newNode = new PHCompositeNode(name.c_str());
  topnodemap[name] = newNode;
  return 0;
}

PHCompositeNode *Fun4AllServer::getNode(const string &name, const string &topnodename)
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
Fun4AllServer::getInputManager(const char *name)
{
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    if ((*iter)->getInputManager(name))
    {
      return (*iter)->getInputManager(name);
    }
  }
  cout << "Could not locate input manager " << name << endl;
  return 0;
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
    cout << "Fun4AllServer: Runnumber forced to " << runnumber << " by RUNNUMBER IntFlag" << endl;
  }
  int iret = 0;
  int icnt = 0;
  int icnt_good = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
  while (!iret)
  {
    int resetnodetree = 0;
    for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
    {
      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << "executing run for input master " << (*iter)->Name() << endl;
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
      if (retval == 1)
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
      ResetNodeTree();
      for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
      {
        (*iter)->PushBackInputMgrsEvents(1);
      }
      continue;
    }
    if (iret)
    {
      break;
    }
    int currentrun = 0;
    for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
    {
      int runno = (*iter)->CurrentRun();
      //	  cout << (*iter)->Name() << " run no: " << runno << endl;
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
            cout << "Mixing of Runs within same event is not supported" << endl;
            cout << "Here is the list of Sync Managers and their runnumbers:" << endl;
            vector<Fun4AllSyncManager *>::const_iterator syiter;
            for (syiter = SyncManagers.begin(); syiter != SyncManagers.end(); ++syiter)
            {
              cout << (*syiter)->Name() << " run number: " << (*syiter)->CurrentRun() << endl;
            }
            cout << "Exiting now" << endl;
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

    if (Verbosity() >= VERBOSITY_SOME)
    {
      // print event cycle counter in log scale if VERBOSITY_SOME
      double significand = 0;
      if (icnt > 0)
      {
        significand = icnt / pow(10, (int) (log10(icnt)));
      }

      if ((fmod(significand, 1.0) == 0 && significand <= 10) or icnt == 0)
      {
        cout << "Fun4AllServer::run - process_event cycle "
             << icnt << "\t for run " << runnumber;
        if (require_nevents)
          cout << ", " << icnt_good << " good event so far";
        cout << endl;
      }
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

    if (require_nevents)
    {
      if (std::find(RetCodes.begin(),
                    RetCodes.end(),
                    static_cast<int>(Fun4AllReturnCodes::ABORTEVENT)) == RetCodes.end())
        icnt_good++;
      if (iret || (nevnts > 0 && icnt_good >= nevnts))
        break;
    }
    else if (iret || (nevnts > 0 && ++icnt >= nevnts))
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
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->skip(nevnts);
  }
  return iret;
}

//_________________________________________________________________
int Fun4AllServer::fileopen(const char *managername, const char *filename)
{
  int iret = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->fileopen(managername, filename);
  }
  return iret;
}

int Fun4AllServer::BranchSelect(const char *managername, const char *branch, int iflag)
{
  int iret = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->BranchSelect(managername, branch, iflag);
  }
  return iret;
}

int Fun4AllServer::BranchSelect(const char *branch, int iflag)
{
  int iret = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->BranchSelect(branch, iflag);
  }
  return iret;
}

int Fun4AllServer::setBranches(const char *managername)
{
  int iret = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->setBranches(managername);
  }
  return iret;
}

int Fun4AllServer::setBranches()
{
  int iret = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    iret += (*iter)->setBranches();
  }
  return iret;
}

int Fun4AllServer::fileclose(const string &managername)
{
  int iret = 0;
  vector<Fun4AllSyncManager *>::const_iterator iter;
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

void Fun4AllServer::GetInputFullFileList(vector<string> &fnames) const
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
  map<string, PHCompositeNode *>::const_iterator iter;
  for (iter = topnodemap.begin(); iter != topnodemap.end(); ++iter)
  {
    names.push_back(iter->first);
  }
  return names.size();
}

void Fun4AllServer::GetOutputManagerList(std::vector<std::string> &names) const
{
  names.clear();
  vector<Fun4AllOutputManager *>::const_iterator iter;
  for (iter = OutputManager.begin(); iter != OutputManager.end(); ++iter)
  {
    names.push_back((*iter)->Name());
  }
  return;
}

void Fun4AllServer::GetModuleList(std::vector<std::string> &names) const
{
  names.clear();
  vector<pair<SubsysReco *, PHCompositeNode *> >::const_iterator iter;
  for (iter = Subsystems.begin(); iter != Subsystems.end(); ++iter)
  {
    names.push_back((*iter).first->Name());
  }
  return;
}

int Fun4AllServer::registerSyncManager(Fun4AllSyncManager *newmaster)
{
  BOOST_FOREACH (Fun4AllSyncManager *syncman, SyncManagers)
  {
    if (syncman->Name() == newmaster->Name())
    {
      cout << "Input Master " << newmaster->Name()
           << " already registered" << endl;
      return -1;
    }
  }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "Registering Input Master " << newmaster->Name() << endl;
  }
  SyncManagers.push_back(newmaster);
  return 0;
}

Fun4AllSyncManager *
Fun4AllServer::getSyncManager(const string &name)
{
  if (name == "DefaultSyncManager")
  {
    return defaultSyncManager;
  }
  vector<Fun4AllSyncManager *>::iterator iter;

  for (iter = SyncManagers.begin(); iter != SyncManagers.end(); ++iter)
  {
    if ((*iter)->Name() == name)
    {
      return *iter;
    }
  }
  cout << "Could not find Input Master " << name << endl;
  return 0;
}

int Fun4AllServer::setRun(const int runno)
{
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runno);
  PHTimeStamp *tstamp = nullptr;
  if (!tstamp)
  {
    tstamp = new PHTimeStamp(0);
    cout << "Fun4AllServer::setRun(): could not get timestamp for run  " << runno
         << ", using tics(0) timestamp: ";
    tstamp->print();
    cout << endl;
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
    cout << "Could not locate node " << name
         << " or no PHObject Node" << endl;
  }
  return;
}

void Fun4AllServer::PrintTimer(const string &name)
{
  map<const string, PHTimer>::const_iterator iter;
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
      cout << "No timer with name " << name << " found" << endl;
      cout << "Existing timers:" << endl;
      for (iter = timer_map.begin(); iter != timer_map.end(); ++iter)
      {
	cout << iter->first << endl;
      }
    }
  }
  return;
}
