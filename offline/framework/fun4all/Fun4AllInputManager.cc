#include "Fun4AllInputManager.h"

#include "Fun4AllServer.h"
#include "SubsysReco.h"

#include <phool/phool.h>

#include <boost/filesystem.hpp>

#include <cstdint>  // for uintmax_t
#include <fstream>
#include <iostream>

Fun4AllInputManager::Fun4AllInputManager(const std::string &name, const std::string &nodename, const std::string &topnodename)
  : Fun4AllBase(name)
  , m_InputNode(nodename)
  , m_TopNodeName(topnodename)
{
  return;
}

Fun4AllInputManager::~Fun4AllInputManager()
{
  while (m_SubsystemsVector.begin() != m_SubsystemsVector.end())
  {
    if (Verbosity())
    {
      m_SubsystemsVector.back()->Verbosity(Verbosity());
    }
    delete m_SubsystemsVector.back();
    m_SubsystemsVector.pop_back();
  }
}

void Fun4AllInputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "SUBSYSTEMS")
  {
    // loop over the map and print out the content (name and location in memory)
    std::cout << "--------------------------------------" << std::endl
              << std::endl;
    std::cout << "List of SubsysRecos in Fun4AllInputManager " << Name() << ":" << std::endl;

    for (SubsysReco *subsys : m_SubsystemsVector)
    {
      std::cout << subsys->Name() << std::endl;
    }
    std::cout << std::endl;
  }
  return;
}

int Fun4AllInputManager::registerSubsystem(SubsysReco *subsystem)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  int iret = subsystem->Init(se->topNode(m_TopNodeName));
  if (iret)
  {
    std::cout << PHWHERE << " Error initializing subsystem "
              << subsystem->Name() << ", return code: " << iret << std::endl;
    return iret;
  }
  if (Verbosity() > 0)
  {
    std::cout << "Registering Subsystem " << subsystem->Name() << std::endl;
  }
  m_SubsystemsVector.push_back(subsystem);
  return 0;
}

int Fun4AllInputManager::RejectEvent()
{
  if (!m_SubsystemsVector.empty())
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    for (SubsysReco *subsys : m_SubsystemsVector)
    {
      if (!m_InitRun)
      {
        subsys->InitRun(se->topNode(m_TopNodeName));
        m_InitRun = 1;
      }
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": Fun4AllInpuManager::EventReject processing " << subsys->Name() << std::endl;
      }
      if (subsys->process_event(se->topNode(m_TopNodeName)) != Fun4AllReturnCodes::EVENT_OK)
      {
        return Fun4AllReturnCodes::DISCARDEVENT;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
