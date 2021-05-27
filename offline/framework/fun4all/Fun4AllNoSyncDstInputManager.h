// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLNOSYNCDSTINPUTMANAGER_H
#define FUN4ALL_FUN4ALLNOSYNCDSTINPUTMANAGER_H

#include "Fun4AllDstInputManager.h"

#include "Fun4AllReturnCodes.h"

#include <string>  // for string

class PHNodeIOManager;
class SyncObject;

class Fun4AllNoSyncDstInputManager : public Fun4AllDstInputManager
{
 public:
  Fun4AllNoSyncDstInputManager(const std::string& name = "DUMMY", const std::string& nodename = "DST", const std::string& topnodename = "TOP");

  ~Fun4AllNoSyncDstInputManager() override {}

  // Effectivly turn off the synchronization checking
  //
  int SyncIt(const SyncObject* /*mastersync*/) override { return Fun4AllReturnCodes::SYNC_OK; }
  int GetSyncObject(SyncObject** /*mastersync*/) override { return Fun4AllReturnCodes::SYNC_NOOBJECT; }
  int NoSyncPushBackEvents(const int nevt) override { return PushBackEvents(nevt); }
  // no sync object we don't need to enable the sync variables
  int setSyncBranches(PHNodeIOManager* /*IManager*/) override { return 0; }

  // turn off reading of the runwise TTree to make run mixing for embedding possible
  int NoRunTTree();

  int SkipForThisManager(const int nevents) override { return PushBackEvents(nevents); }
  int HasSyncObject() const override { return 0; }
};

#endif  // FUN4ALL_FUN4ALLNOSYNCDSTINPUTMANAGER_H
