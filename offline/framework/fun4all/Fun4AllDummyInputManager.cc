#include "Fun4AllDummyInputManager.h"

#include "Fun4AllSyncManager.h"

#include <phool/recoConsts.h>

#include <iostream>

Fun4AllDummyInputManager::Fun4AllDummyInputManager(const std::string &name, const std::string &nodename)
  : Fun4AllInputManager(name, nodename)
{
  FileName("NOFILE-0000000000-0000.root");
  return;
}

int Fun4AllDummyInputManager::ResetFileList()
{
  m_NumEvents = 0;
  return 0;
}

int Fun4AllDummyInputManager::PushBackEvents(const int nevt)
{
  m_NumEvents -= nevt;
  m_SumEvents -= nevt;
  return 0;
}

void Fun4AllDummyInputManager::setSyncManager(Fun4AllSyncManager *master)
{
  // set the runnumber in Fun4All if rc flag is set
  // so InitRun is triggered.
  // This setSyncManager is executed during Fun4AllServer::registerInputManager()
  // normally the runnumber is set in Fun4AllInputManager::fileopen() but since it
  // would be kind of ridicolous to call this for a dummy input manager
  // we set the runnumber here
  Fun4AllInputManager::setSyncManager(master);
  recoConsts *rc = recoConsts::instance();
  int runnumber = rc->get_IntFlag("RUNNUMBER", 0);
  MySyncManager()->CurrentRun(runnumber);
  return;
}

int Fun4AllDummyInputManager::run(const int nevents)
{
  m_NumEvents += nevents;
  m_SumEvents += nevents;
  MySyncManager()->CurrentEvent(m_NumEvents);
  if (Verbosity() > 0)
  {
    std::cout << "Event No: " << m_NumEvents;
    if (m_SumEvents != m_NumEvents)
    {
      std::cout << ", Event Sum: " << m_SumEvents;
    }
    std::cout << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
