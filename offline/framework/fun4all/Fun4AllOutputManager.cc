#include "Fun4AllOutputManager.h"

#include "Fun4AllServer.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

Fun4AllOutputManager::Fun4AllOutputManager(const string &name)
  : Fun4AllBase(name)
  , m_NEvents(0)
{
}

Fun4AllOutputManager::Fun4AllOutputManager(const string &name, const string &outfname)
  : Fun4AllBase(name)
  , m_NEvents(0)
  , m_OutFileName(outfname)
{
}

//___________________________________________________________________
int Fun4AllOutputManager::AddEventSelector(const string &recomodule)
{
  string newselector = recomodule;
  for (string evtsel : m_EventSelectorsVector)
  {
    if (evtsel == newselector)
    {
      cout << "Event Selector " << newselector << " allready in list" << endl;
      return -1;
    }
  }
  cout << "EventSelector: " << &m_EventSelectorsVector << endl;
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
void Fun4AllOutputManager::Print(const string &what) const
{
  if (what == "ALL" || what == "EVENTSELECTOR")
  {
    unsigned icnt = 0;
    for (string evtsel : m_EventSelectorsVector)
    {
      cout << Name() << ": Reco Module " << evtsel << " select Events" << endl;
      cout << Name() << ": Reco Module Index: " << m_RecoModuleIndexVector[icnt] << endl;
      icnt++;
    }
  }
  if (what == "ALL" || what == "EVENTSWRITTEN")
  {
    cout << Name() << " wrote " << EventsWritten() << " Events" << endl;
  }
  return;
}

//___________________________________________________________________
int Fun4AllOutputManager::DoNotWriteEvent(vector<int> *retcodes) const
{
  int iret = 0;
  for (unsigned int index : m_RecoModuleIndexVector)
  {
    iret += (*retcodes)[index];
  }
  return iret;
}
