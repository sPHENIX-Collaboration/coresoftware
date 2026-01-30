// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLINPUTMANAGER_H
#define FUN4ALL_FUN4ALLINPUTMANAGER_H

#include "Fun4AllBase.h"
#include "Fun4AllReturnCodes.h"
#include "InputFileHandler.h"

#include <list>
#include <string>
#include <type_traits>  // for __decay_and_strip<>::__type
#include <utility>      // for make_pair, pair
#include <vector>

class PHCompositeNode;
class SubsysReco;
class SyncObject;
class Fun4AllSyncManager;

class Fun4AllInputManager : public Fun4AllBase, public InputFileHandler
{
 public:
  using Fun4AllBase::Verbosity;
  ~Fun4AllInputManager() override;
  virtual int run(const int /*nevents=0*/) { return -1; }
  virtual int ReadInRunNode(PHCompositeNode * /*RunNode*/) { return -1; }
  virtual int GetSyncObject(SyncObject ** /*mastersync*/) { return 0; }
  virtual int SyncIt(const SyncObject * /*mastersync*/) { return Fun4AllReturnCodes::SYNC_FAIL; }
  virtual int BranchSelect(const std::string & /*branch*/, const int /*iflag*/) { return -1; }
  virtual int setBranches() { return -1; }  // publich bc needed by the sync manager
  virtual int SkipForThisManager(const int /*nevents*/) { return 0; }
  virtual int HasSyncObject() const { return 0; }
  virtual std::string GetString(const std::string &) const { return ""; }
  virtual int PushBackEvents(const int /*nevt*/) { return -1; }
  virtual int RejectEvent();
  // so people can use the skip they are used to instead of PushBackEvents
  // with negative arg
  virtual int skip(const int nevt) { return PushBackEvents(-nevt); }
  virtual int NoSyncPushBackEvents(const int /*nevt*/) { return -1; }
  virtual void setSyncManager(Fun4AllSyncManager *master) { m_MySyncManager = master; }
  virtual int ResetEvent() { return 0; }
  virtual void SetRunNumber(const int runno) { m_MyRunNumber = runno; }
  virtual int RunNumber() const { return m_MyRunNumber; }

  void Print(const std::string &what = "ALL") const override;

  int registerSubsystem(SubsysReco *subsystem);
  const std::string &InputNode() { return m_InputNode; }
  void InputNode(const std::string &innode) { m_InputNode = innode; }
  const std::string &TopNodeName() const { return m_TopNodeName; }
  void Verbosity(const uint64_t ival) override;
  
 protected:
  Fun4AllInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  Fun4AllSyncManager *MySyncManager() { return m_MySyncManager; }
  void DisableReadCache() { m_disable_read_cache_flag = true; }
  bool ReadCacheDisabled() const { return m_disable_read_cache_flag; }

 private:
  Fun4AllSyncManager *m_MySyncManager {nullptr};
  int m_MyRunNumber {0};
  int m_InitRun {0};
  bool m_disable_read_cache_flag{false};
  std::vector<SubsysReco *> m_SubsystemsVector;
  std::string m_InputNode;
  std::string m_TopNodeName;
};

#endif
