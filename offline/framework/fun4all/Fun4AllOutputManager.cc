#include "Fun4AllOutputManager.h"

#include <iostream>
#include <string>
#include <vector>

Fun4AllOutputManager::Fun4AllOutputManager(const std::string &name)
  : Fun4AllBase(name)
{
}

Fun4AllOutputManager::Fun4AllOutputManager(const std::string &name, const std::string &outfname)
  : Fun4AllBase(name)
  , m_OutFileName(outfname)
{
}

//___________________________________________________________________
int Fun4AllOutputManager::AddEventSelector(const std::string &recomodule)
{
  const std::string &newselector = recomodule;
  for (const std::string &evtsel : m_EventSelectorsVector)
  {
    if (evtsel == newselector)
    {
      std::cout << "Event Selector " << newselector << " allready in list" << std::endl;
      return -1;
    }
  }
  std::cout << "EventSelector: " << &m_EventSelectorsVector << std::endl;
  m_EventSelectorsVector.push_back(newselector);
  return 0;
}

//___________________________________________________________________
int Fun4AllOutputManager::WriteGeneric(PHCompositeNode *startNode)
{
  m_NEvents++;
  int iret = Write(startNode);
  return iret;
}

//___________________________________________________________________
void Fun4AllOutputManager::Print(const std::string &what) const
{
  if (what == "ALL" || what == "EVENTSELECTOR")
  {
    unsigned icnt = 0;
    for (const std::string &evtsel : m_EventSelectorsVector)
    {
      std::cout << Name() << ": Reco Module " << evtsel << " select Events" << std::endl;
      std::cout << Name() << ": Reco Module Index: " << m_RecoModuleIndexVector[icnt] << std::endl;
      icnt++;
    }
  }
  if (what == "ALL" || what == "EVENTSWRITTEN")
  {
    std::cout << Name() << " wrote " << EventsWritten() << " Events" << std::endl;
  }
  return;
}

//___________________________________________________________________
int Fun4AllOutputManager::DoNotWriteEvent(std::vector<int> *retcodes) const
{
  int iret = 0;
  for (unsigned int index : m_RecoModuleIndexVector)
  {
    iret += (*retcodes)[index];
  }
  return iret;
}
