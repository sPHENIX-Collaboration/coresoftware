#ifndef FUN4ALLDSTINPUTMANAGER_H__
#define FUN4ALLDSTINPUTMANAGER_H__

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
  int isOpen() { return isopen; }
  int GetSyncObject(SyncObject **mastersync);
  int SyncIt(const SyncObject *mastersync);
  int BranchSelect(const std::string &branch, const int iflag);
  int setBranches();
  virtual int setSyncBranches(PHNodeIOManager *IManager);
  void Print(const std::string &what = "ALL") const;
  int PushBackEvents(const int i);

 protected:
  int ReadNextEventSyncObject();
  int OpenNextFile();
  int readrunttree;
  int isopen;
  int events_total;
  int events_thisfile;
  int events_skipped_during_sync;
  std::string fullfilename;
  std::string RunNode;
  std::map<const std::string, int> branchread;
  std::string syncbranchname;
  PHCompositeNode *dstNode;
  PHCompositeNode *runNode;
  PHCompositeNode *runNodeCopy;
  PHCompositeNode *runNodeSum;
  PHNodeIOManager *IManager;
  SyncObject *syncobject;
};

#endif /* __FUN4ALLDSTINPUTMANAGER_H__ */
