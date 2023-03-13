// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLINPUTMANAGER_H
#define FUN4ALL_FUN4ALLINPUTMANAGER_H

#include "Fun4AllBase.h"
#include "Fun4AllReturnCodes.h"

#include <list>
#include <string>
#include <type_traits>  // for __decay_and_strip<>::__type
#include <utility>      // for make_pair, pair
#include <vector>

class PHCompositeNode;
class SubsysReco;
class SyncObject;
class Fun4AllSyncManager;

class Fun4AllInputManager : public Fun4AllBase
{
 public:
  ~Fun4AllInputManager() override;
  virtual int fileopen(const std::string & /*filename*/) { return -1; }
  virtual int fileclose() { return -1; }
  virtual int run(const int /*nevents=0*/) { return -1; }
  virtual int ReadInRunNode(PHCompositeNode * /*RunNode*/) { return -1; }
  std::string FileName() const { return m_FileName; }
  void FileName(const std::string &fn) { m_FileName = fn; }
  virtual int GetSyncObject(SyncObject ** /*mastersync*/) { return 0; }
  virtual int SyncIt(const SyncObject * /*mastersync*/) { return Fun4AllReturnCodes::SYNC_FAIL; }
  virtual int BranchSelect(const std::string & /*branch*/, const int /*iflag*/) { return -1; }
  virtual int setBranches() { return -1; }  // publich bc needed by the sync manager
  void Print(const std::string &what = "ALL") const override;
  virtual int PushBackEvents(const int /*nevt*/) { return -1; }
  // so people can use the skip they are used to instead of PushBackEvents
  // with negative arg
  virtual int skip(const int nevt) { return PushBackEvents(-nevt); }
  virtual int NoSyncPushBackEvents(const int /*nevt*/) { return -1; }
  int AddFile(const std::string &filename);
  int AddListFile(const std::string &filename, const int do_it = 0);
  int registerSubsystem(SubsysReco *subsystem);
  virtual int RejectEvent();
  void Repeat(const int i = -1) { m_Repeat = i; }
  virtual void setSyncManager(Fun4AllSyncManager *master) { m_MySyncManager = master; }
  virtual int ResetFileList();
  virtual int ResetEvent() { return 0; }
  virtual void SetRunNumber(const int runno) { m_MyRunNumber = runno; }
  virtual int RunNumber() const { return m_MyRunNumber; }
  void AddToFileOpened(const std::string &filename) { m_FileListOpened.push_back(filename); }
  std::pair<std::list<std::string>::const_iterator, std::list<std::string>::const_iterator> FileOpenListBeginEnd() { return std::make_pair(m_FileListOpened.begin(), m_FileListOpened.end()); }
  std::string InputNode() { return m_InputNode; }
  void InputNode(const std::string &innode) { m_InputNode = innode; }
  std::string TopNodeName() const { return m_TopNodeName; }
  bool FileListEmpty() const { return m_FileList.empty(); }
  virtual int IsOpen() const { return m_IsOpen; }
  virtual int SkipForThisManager(const int /*nevents*/) { return 0; }
  virtual int HasSyncObject() const { return 0; }
  virtual std::string GetString(const std::string &) const { return ""; }

 protected:
  Fun4AllInputManager(const std::string &name = "DUMMY", const std::string &nodename = "DST", const std::string &topnodename = "TOP");
  void UpdateFileList();
  int OpenNextFile();
  void IsOpen(const int i) { m_IsOpen = i; }
  Fun4AllSyncManager *MySyncManager() { return m_MySyncManager; }

 private:
  Fun4AllSyncManager *m_MySyncManager = nullptr;
  int m_IsOpen = 0;
  int m_Repeat = 0;
  int m_MyRunNumber = 0;
  int m_InitRun = 0;
  std::vector<SubsysReco *> m_SubsystemsVector;
  std::string m_InputNode;
  std::string m_FileName;
  std::string m_TopNodeName;
  std::list<std::string> m_FileList;
  std::list<std::string> m_FileListCopy;
  std::list<std::string> m_FileListOpened;  // all files which were opened during running
};

#endif
