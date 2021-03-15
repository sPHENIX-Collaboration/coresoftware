// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLDSTINPUTMANAGER_H
#define FUN4ALL_FUN4ALLDSTINPUTMANAGER_H

#include "Fun4AllInputManager.h"

#include <map>
#include <string>

class PHCompositeNode;
class PHNodeIOManager;
class SyncObject;

class Fun4AllDstInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllDstInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  virtual ~Fun4AllDstInputManager();
  int fileopen(const std::string &filenam);
  int fileclose();
  int run(const int nevents = 0);
  int GetSyncObject(SyncObject **mastersync);
  int SyncIt(const SyncObject *mastersync);
  int BranchSelect(const std::string &branch, const int iflag);
  int setBranches();
  virtual int setSyncBranches(PHNodeIOManager *IManager);
  void Print(const std::string &what = "ALL") const;
  int PushBackEvents(const int i);
  virtual int HasSyncObject() const;

 protected:
  int ReadNextEventSyncObject();
  void ReadRunTTree(const int i) { m_ReadRunTTree = i; }

 private:
  int m_ReadRunTTree = 1;
  int events_total = 0;
  int events_thisfile = 0;
  int events_skipped_during_sync = 0;
  int m_HaveSyncObject = 0;
  std::string fullfilename;
  std::string RunNode = "RUN";
  std::map<const std::string, int> branchread;
  std::string syncbranchname;
  PHCompositeNode *dstNode = nullptr;
  PHCompositeNode *runNode = nullptr;
  PHCompositeNode *runNodeCopy = nullptr;
  PHCompositeNode *runNodeSum = nullptr;
  PHNodeIOManager *IManager = nullptr;
  SyncObject *syncobject = nullptr;
};

#endif /* __FUN4ALLDSTINPUTMANAGER_H__ */
