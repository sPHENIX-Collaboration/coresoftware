#include "Fun4AllDummyInputManager.h"
#include "Fun4AllSyncManager.h"

#include <phool/recoConsts.h>

#include <iostream>

using namespace std;

Fun4AllDummyInputManager::Fun4AllDummyInputManager(const string &name, const string &nodename): 
  Fun4AllInputManager(name, nodename),
  numevents(0)
{
  FileName("NOFILE-0000000000-0000.root");
  return;
}

void
Fun4AllDummyInputManager::setSyncManager(Fun4AllSyncManager *master)
{
  // set the runnumber in Fun4All if rc flag is set
  // so InitRun is triggered. 
  // This setSyncManager is executed during Fun4AllServer::registerInputManager()
  // normally the runnumber is set in Fun4AllInputManager::fileopen() but since it
  // would be kind of ridicolous to call this for a dummy input manager
  // we set the runnumber here
  Fun4AllInputManager::setSyncManager(master);
  recoConsts *rc = recoConsts::instance();
  int runnumber = rc->get_IntFlag("RUNNUMBER",0);
  mySyncManager->CurrentRun(runnumber);
  return;
}

int
Fun4AllDummyInputManager::run(const int nevents)
{
  numevents+= nevents;
  if (Verbosity()>0)
    {
      cout << "Event No: " << numevents << endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
