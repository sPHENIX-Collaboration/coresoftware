// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLDUMMYINPUTMANAGER_H
#define FUN4ALL_FUN4ALLDUMMYINPUTMANAGER_H

#include "Fun4AllInputManager.h"

#include "Fun4AllReturnCodes.h"

#include <string>  // for string

class Fun4AllSyncManager;
class SyncObject;

class Fun4AllDummyInputManager : public Fun4AllInputManager
{
 public:
  Fun4AllDummyInputManager(const std::string& name = "DUMMY", const std::string& nodename = "DST");
  virtual ~Fun4AllDummyInputManager() {}
  int fileopen(const std::string&) { return 0; }
  int fileclose() { return 0; }
  int IsOpen() const { return 1; }
  int run(const int /*nevents=0*/);
  int GetSyncObject(SyncObject** /*mastersync*/) { return Fun4AllReturnCodes::SYNC_NOOBJECT; }
  int SyncIt(const SyncObject* /*mastersync*/) { return Fun4AllReturnCodes::SYNC_OK; }
  void setSyncManager(Fun4AllSyncManager* master);
  int PushBackEvents(const int nevt);
  int NoSyncPushBackEvents(const int nevt) { return PushBackEvents(nevt); }
  int ResetFileList();

 private:
  int m_NumEvents = 0;
  int m_SumEvents = 0;
};

#endif
