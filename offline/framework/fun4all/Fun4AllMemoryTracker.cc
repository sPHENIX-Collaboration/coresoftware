#include "Fun4AllMemoryTracker.h"

#include <TSystem.h>

#include <iostream>

using namespace std;

Fun4AllMemoryTracker *Fun4AllMemoryTracker::mInstance = nullptr;

Fun4AllMemoryTracker::Fun4AllMemoryTracker()
  : Fun4AllBase("Fun4AllMemoryTracker")
{
}

int Fun4AllMemoryTracker::GetRSSMemory() const
{
  ProcInfo_t procinfo;
  gSystem->GetProcInfo(&procinfo);
  return procinfo.fMemResident;
}

Fun4AllMemoryTracker::~Fun4AllMemoryTracker()
{
  mMemoryTrackerMap.clear();
  mStartMem.clear();
}

void Fun4AllMemoryTracker::Snapshot(const string &trackername, const string &group)
{
  string name = CreateFullTrackerName(trackername, group);
  auto iter = mMemoryTrackerMap.find(name);
  if (iter != mMemoryTrackerMap.end())
  {
    iter->second.push_back(GetRSSMemory());
  }
  else
  {
    vector<int> mvec;
    mvec.push_back(GetRSSMemory());
    mMemoryTrackerMap.insert(make_pair(name, mvec));
  }
  if (Verbosity() > 0)
  {
    cout << "Snapshot name: " << name << ", mem: " << GetRSSMemory() << endl;
  }
  return;
}

void Fun4AllMemoryTracker::Start(const string &trackername, const string &group)
{
  string name = CreateFullTrackerName(trackername, group);
  auto iter = mStartMem.find(name);
  int RSSMemory = GetRSSMemory();
  if (iter != mStartMem.end())
  {
    iter->second = RSSMemory;
  }
  else
  {
    mStartMem.insert(make_pair(name, RSSMemory));
  }
  if (Verbosity() > 0)
  {
    cout << "Start name: " << name << ", mem: " << RSSMemory << endl;
  }
}

void Fun4AllMemoryTracker::Stop(const string &trackername, const string &group)
{
  string name = CreateFullTrackerName(trackername, group);
  auto iter = mStartMem.find(name);
  int RSSMemory = GetRSSMemory();
  if (iter != mStartMem.end())
  {
    int diff = RSSMemory - iter->second;
    auto iterM = mMemoryTrackerMap.find(name);
    if (iterM != mMemoryTrackerMap.end())
    {
      iterM->second.push_back(diff);
    }
    else
    {
      vector<int> mvec;
      mvec.push_back(diff);
      mMemoryTrackerMap.insert(make_pair(name, mvec));
    }
    if (Verbosity() > 0)
    {
      cout << "Stop name: " << name << ", mem: " << RSSMemory << ", diff: " << diff << endl;
    }
  }
  return;
}

string Fun4AllMemoryTracker::CreateFullTrackerName(const string &trackername, const string &group)
{
  string name = trackername;
  if (!group.empty())
  {
    name = group + "_" + name;
  }
  return name;
}

void Fun4AllMemoryTracker::PrintMemoryTracker(const string &name) const
{
  map<const string, std::vector<int>>::const_iterator iter;
  if (name.empty())
  {
    for (iter = mMemoryTrackerMap.begin(); iter != mMemoryTrackerMap.end(); ++iter)
    {
      cout << "SubsysReco/OutputManager: " << iter->first << endl;
      vector<int> memvec = iter->second;
      for (auto vit = memvec.begin(); vit != memvec.end(); ++vit)
      {
        cout << *vit << " ";
      }
      cout << endl;
    }
  }
  else
  {
    iter = mMemoryTrackerMap.find(name);
    if (iter != mMemoryTrackerMap.end())
    {
      cout << "SubsysReco/OutputManager: " << iter->first << endl;
      vector<int> memvec = iter->second;
      for (auto vit = memvec.begin(); vit != memvec.end(); ++vit)
      {
        cout << *vit << " ";
      }
      cout << endl;
    }
    else
    {
      cout << "No Memory Tracker with name " << name << " found" << endl;
      cout << "Existing Memory Trackers:" << endl;
      for (iter = mMemoryTrackerMap.begin(); iter != mMemoryTrackerMap.end(); ++iter)
      {
        cout << iter->first << endl;
      }
    }
  }
  return;
}

vector<int> Fun4AllMemoryTracker::GetMemoryVector(const std::string &name) const
{
  vector<int> memvec;
  auto iter =  mMemoryTrackerMap.find(name);
  if (iter != mMemoryTrackerMap.end())
  {
    memvec = iter->second;
  }
  return memvec;
}
