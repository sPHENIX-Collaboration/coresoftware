#ifndef FUN4ALLSERVER_H
#define FUN4ALLSERVER_H

#include "Fun4AllBase.h"
#include "Fun4AllHistoManager.h"

#include <phool/PHTimer.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

class Fun4AllInputManager;
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
  virtual ~Fun4AllServer();

  virtual bool registerHisto(const char *hname, TNamed *h1d, const int replace = 0);
  virtual bool registerHisto(TNamed *h1d, const int replace = 0);
  template <typename T>
  T *makeHisto(T *t)
  {
    return ServerHistoManager->makeHisto(t);
  }
  virtual int isHistoRegistered(const std::string &name) const;

  int registerSubsystem(SubsysReco *subsystem, const std::string &topnodename = "TOP");
  int unregisterSubsystem(SubsysReco *subsystem);
  SubsysReco *getSubsysReco(const std::string &name);
  int registerOutputManager(Fun4AllOutputManager *manager);
  Fun4AllOutputManager *getOutputManager(const std::string &name);
  int registerHistoManager(Fun4AllHistoManager *manager);
  Fun4AllHistoManager *getHistoManager(const std::string &name);
  TNamed *getHisto(const std::string &hname) const;
  TNamed *getHisto(const unsigned int ihisto) const;
  const char *getHistoName(const unsigned int ihisto) const;
  virtual void Print(const std::string &what = "ALL") const;

  void InitAll();
  int BeginRunTimeStamp(PHTimeStamp &TimeStp);
  int dumpHistos(const std::string &filename = "", const std::string &openmode = "RECREATE");
  int process_event(PHCompositeNode *topNode);
  int Reset();
  virtual int BeginRun(const int runno);
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
  Fun4AllInputManager *getInputManager(const char *name);
  int PrdfEvents() const;
  int DstEvents() const;

  //! run n events (0 means up to end of file)
  int run(const int nevnts = 0, const bool require_nevents = false);

  /*! 
    \brief skip n events (0 means up to the end of file). 
    Skip means read, don't process.
  */
  int skip(const int nevnts = 0);

  int fileopen(const char *managername = "NONE", const char *filename = "NONE");
  int fileclose(const std::string &managername = "");
  int SegmentNumber();
  int ResetNodeTree();
  int BranchSelect(const char *managername, const char *branch, int iflag);
  int BranchSelect(const char *branch, int iflag);
  int setBranches(const char *managername);
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
  void NodeIdentify(const std::string &name);
  void KeepDBConnection(const int i = 1) { keep_db_connected = i; }
  void PrintTimer(const std::string &name = "");

 protected:
  Fun4AllServer(const std::string &name = "Fun4AllServer");
  int InitNodeTree(PHCompositeNode *topNode);
  int CountOutNodes(PHCompositeNode *startNode);
  int CountOutNodesRecursive(PHCompositeNode *startNode, const int icount);
  int UpdateEventSelector(Fun4AllOutputManager *manager);
  int unregisterSubsystemsNow();
  int setRun(const int runnumber);
  static Fun4AllServer *__instance;
  int OutNodeCount;
  int bortime_override;
  int ScreamEveryEvent;
  int unregistersubsystem;
  int runnumber;
  int eventnumber;
  std::vector<std::string> ComplaintList;
  PHCompositeNode *TopNode;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *> > Subsystems;
  std::vector<std::pair<SubsysReco *, PHCompositeNode *> > DeleteSubsystems;
  std::vector<int> RetCodes;
  std::vector<Fun4AllOutputManager *> OutputManager;
  std::vector<TDirectory *> TDirCollection;
  Fun4AllHistoManager *ServerHistoManager;
  std::vector<Fun4AllHistoManager *> HistoManager;
  std::map<std::string, PHCompositeNode *> topnodemap;
  PHTimeStamp *beginruntimestamp;
  std::string default_Tdirectory;
  Fun4AllSyncManager *defaultSyncManager;
  std::vector<Fun4AllSyncManager *> SyncManagers;
  std::map<int, int> retcodesmap;
  std::map<const std::string, PHTimer> timer_map;
  TH1 *FrameWorkVars;
  int keep_db_connected;
};

#endif /* __FUN4ALLSERVER_H */
