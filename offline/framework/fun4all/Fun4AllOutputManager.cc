#include "Fun4AllOutputManager.h"
#include "Fun4AllServer.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

Fun4AllOutputManager::Fun4AllOutputManager(const string &name)
  : Fun4AllBase(name)
  , nEvents(0)
{
}

//___________________________________________________________________
int Fun4AllOutputManager::AddEventSelector(const string &recomodule)
{
  string newselector = recomodule;
  vector<string>::iterator iter;

  for (iter = EventSelectors.begin(); iter != EventSelectors.end(); ++iter)
    if (*iter == newselector)
    {
      cout << "Event Selector " << newselector << " allready in list" << endl;
      return -1;
    }

  cout << "EventSelector: " << &EventSelectors << endl;
  EventSelectors.push_back(newselector);
  return 0;
}

//___________________________________________________________________
int Fun4AllOutputManager::WriteGeneric(PHCompositeNode *startNode)
{
  nEvents++;
  int iret = Write(startNode);
  return iret;
}

//___________________________________________________________________
void Fun4AllOutputManager::Print(const string &what) const
{
  if (what == "ALL" || what == "EVENTSELECTOR")
  {
    unsigned icnt = 0;
    for (vector<string>::const_iterator iter = EventSelectors.begin(); iter != EventSelectors.end(); ++iter)
    {
      cout << Name() << ": Reco Module " << *iter << " select Events" << endl;
      cout << Name() << ": Reco Module Index: " << recomoduleindex[icnt] << endl;
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
  for (vector<unsigned>::const_iterator iter = recomoduleindex.begin(); iter != recomoduleindex.end(); ++iter)
  {
    const unsigned index = *iter;
    iret += (*retcodes)[index];
  }
  return iret;
}
