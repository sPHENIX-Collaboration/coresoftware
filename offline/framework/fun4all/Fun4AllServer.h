// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLSERVER_H
#define FUN4ALL_FUN4ALLSERVER_H

#include "Fun4AllBase.h"

#include "Fun4AllHistoManager.h"  // for Fun4AllHistoManager

#include <phool/PHTimer.h>

#include <deque>
#include <iostream>
#include <map>
#include <string>
#include <utility>  // for pair
#include <vector>

class Fun4AllInputManager;
class Fun4AllMemoryTracker;
class Fun4AllSyncManager;
class Fun4AllOutputManager;
class PHCompositeNode;
class PHTimeStamp;
class SubsysReco;
class TDirectory;
class TH1;
class TNamed;

class Fun4AllServer : public Fun4AllBase
{
 public:
  static Fun4AllServer *instance();
  ~Fun4AllServer() override;

  virtual bool registerHisto(const std::string &hname, TNamed *h1d, const int replace = 0);
  virtual bool registerHisto(TNamed *h1d, const int replace = 0);
  template <typename T>
  T *makeHisto(T *t)
  {
    return ServerHistoManager->makeHisto(t);
  }
  virtual int isHistoRegistered(const std::string &name) const;

  int registerSubsystem(SubsysReco *subsystem, const std::string &topnodename = "TOP");
  void addNewSubsystem(SubsysReco *subsystem, const std::string &topnodename = "TOP") { NewSubsystems.push_back(std::make_pair(subsystem, topnodename)); }
  int unregisterSubsystem(SubsysReco *subsystem);
  SubsysReco *getSubsysReco(const std::string &name);
  int registerOutputManager(Fun4AllOutputManager *manager);
  Fun4AllOutputManager *getOutputManager(const std::string &name);
  int registerHistoManager(Fun4AllHistoManager *manager);
  Fun4AllHistoManager *getHistoManager(const std::string &name);
  TNamed *getHisto(const std::string &hname) const;
  TNamed *getHisto(const unsigned int ihisto) const;
  std::string getHistoName(const unsigned int ihisto) const;
  void Print(const std::string &what = "ALL") const override;

  void InitAll();
  int BeginRunTimeStamp(PHTimeStamp &TimeStp);
  int dumpHistos(const std::string &filename, const std::string &openmode = "RECREATE");
  int Reset();
  virtual int BeginRun(const int runno);
  int BeginRunSubsystem(const std::pair<SubsysReco *, PHCompositeNode *> &subsys);
  virtual int EndRun(const int runno = 0);
  virtual int End();
  PHCompositeNode *topNode() const { return TopNode; }
  PHCompositeNode *topNode(const std::string &name);
  int outfileclose();
  virtual int process_event();
  PHCompositeNode *getNode(const std::string &name, const std::string &topnodename = "TOP");
  int AddTopNode(const std::string &name);
  int MakeNodesTransient(PHCompositeNode *startNode);
  int MakeNodesPersistent(PHCompositeNode *startNode);

  int AddComplaint(const std::string &complaint, const std::string &remedy);

  // Interface to the default Input Master
  int registerInputManager(Fun4AllInputManager *InManager);
  Fun4AllInputManager *getInputManager(const std::string &name);
  int PrdfEvents() const;
  int DstEvents() const;

  //! run n events (0 means up to end of file)
  int run(const int nevnts = 0, const bool require_nevents = false);

  /*! 
    \brief skip n events (0 means up to the end of file). 
    Skip means read, don't process.
  */
  int skip(const int nevnts = 0);

  int fileopen(const std::string &managername, const std::string &filename);
  int fileclose(const std::string &managername);
  int SegmentNumber();
  int ResetNodeTree();
  int BranchSelect(const std::string &managername, const std::string &branch, int iflag);
  int BranchSelect(const std::string &branch, int iflag);
  int setBranches(const std::string &managername);
  int setBranches();
  virtual int DisconnectDB();
  virtual void identify(std::ostream &out = std::cout) const;
  unsigned GetTopNodes(std::vector<std::string> &names) const;
  void GetInputFullFileList(std::vector<std::string> &fnames) const;
  void GetOutputManagerList(std::vector<std::string> &names) const;
  void GetModuleList(std::vector<std::string> &names) const;
  Fun4AllSyncManager *getSyncManager(const std::string &name = "DefaultSyncManager");
  int registerSyncManager(Fun4AllSyncManager *newmaster);
  int retcodestats(const int iret) { return retcodesmap[iret]; }
  void EventNumber(const int evtno) { eventnumber = evtno; }
  int EventNumber() const { return eventnumber; }
  void NodeIdentify(const std::string &name);
  void KeepDBConnection(const int i = 1) { keep_db_connected = i; }
  void PrintTimer(const std::string &name = "");
  void PrintMemoryTracker(const std::string &name = "") const;
  int RunNumber() const { return runnumber; }
  int EventCounter() const { return eventcounter; }

 protected:
  Fun4AllServer(const std::string &name = "Fun4AllServer");
  int InitNodeTree(PHCompositeNode *topNode);
  int CountOutNodes(PHCompositeNode *startNode);
  int CountOutNodesRecursive(PHCompositeNode *startNode, const int icount);
  int UpdateEventSelector(Fun4AllOutputManager *manager);
  int unregisterSubsystemsNow();
  int setRun(const int runnumber);
  static Fun4AllServer *__instance;
  TH1 *FrameWorkVars = nullptr;
  Fun4AllMemoryTracker *ffamemtracker = nullptr;
  Fun4AllHistoManager *ServerHistoManager = nullptr;
  PHTimeStamp *beginruntimestamp = nullptr;
  PHCompositeNode *TopNode = nullptr;
  Fun4AllSyncManager *defaultSyncManager = nullptr;

  int OutNodeCount = 0;
  int bortime_override = 0;
  int ScreamEveryEvent = 0;
  int unregistersubsystem = 0;
  int runnumber = 0;
  int eventnumber = 0;
  int eventcounter = 0;
  int keep_db_connected = 0;

  std::vector<std::string> ComplaintList;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>> Subsystems;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *>> DeleteSubsystems;
  std::deque<std::pair<SubsysReco *, std::string>> NewSubsystems;
  std::vector<int> RetCodes;
  std::vector<Fun4AllOutputManager *> OutputManager;
  std::vector<TDirectory *> TDirCollection;
  std::vector<Fun4AllHistoManager *> HistoManager;
  std::map<std::string, PHCompositeNode *> topnodemap;
  std::string default_Tdirectory;
  std::vector<Fun4AllSyncManager *> SyncManagers;
  std::map<int, int> retcodesmap;
  std::map<const std::string, PHTimer> timer_map;
};

#endif
