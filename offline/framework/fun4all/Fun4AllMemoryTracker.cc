#include "Fun4AllMemoryTracker.h"

#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair, make_pair

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

void Fun4AllMemoryTracker::Snapshot(const std::string &trackername, const std::string &group)
{
  std::string name = CreateFullTrackerName(trackername, group);
  auto iter = mMemoryTrackerMap.find(name);
  if (iter != mMemoryTrackerMap.end())
  {
    iter->second.push_back(GetRSSMemory());
  }
  else
  {
    std::vector<int> mvec;
    mvec.push_back(GetRSSMemory());
    mMemoryTrackerMap.insert(make_pair(name, mvec));
  }
  if (Verbosity() > 0)
  {
    std::cout << "Snapshot name: " << name << ", mem: " << GetRSSMemory() << std::endl;
  }
  return;
}

void Fun4AllMemoryTracker::Start(const std::string &trackername, const std::string &group)
{
  std::string name = CreateFullTrackerName(trackername, group);
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
    std::cout << "Start name: " << name << ", mem: " << RSSMemory << std::endl;
  }
}

void Fun4AllMemoryTracker::Stop(const std::string &trackername, const std::string &group)
{
  std::string name = CreateFullTrackerName(trackername, group);
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
      std::vector<int> mvec;
      mvec.push_back(diff);
      mMemoryTrackerMap.insert(make_pair(name, mvec));
    }
    if (Verbosity() > 0)
    {
      std::cout << "Stop name: " << name << ", mem: " << RSSMemory << ", diff: " << diff << std::endl;
    }
  }
  return;
}

std::string Fun4AllMemoryTracker::CreateFullTrackerName(const std::string &trackername, const std::string &group)
{
  std::string name = trackername;
  if (!group.empty())
  {
    name = group + "_" + name;
  }
  return name;
}

void Fun4AllMemoryTracker::PrintMemoryTracker(const std::string &name) const
{
  std::map<std::string, std::vector<int>>::const_iterator iter;
  if (name.empty())
  {
    for (iter = mMemoryTrackerMap.begin(); iter != mMemoryTrackerMap.end(); ++iter)
    {
      std::cout << iter->first << ": ";
      std::vector<int> memvec = iter->second;
      for (int & vit : memvec)
      {
        std::cout << vit << " ";
      }
      std::cout << std::endl;
    }
  }
  else
  {
    iter = mMemoryTrackerMap.find(name);
    if (iter != mMemoryTrackerMap.end())
    {
      std::cout << "SubsysReco/OutputManager: " << iter->first << std::endl;
      std::vector<int> memvec = iter->second;
      for (int & vit : memvec)
      {
        std::cout << vit << " ";
      }
      std::cout << std::endl;
    }
    else
    {
      std::cout << "No Memory Tracker with name " << name << " found" << std::endl;
      std::cout << "Existing Memory Trackers:" << std::endl;
      for (iter = mMemoryTrackerMap.begin(); iter != mMemoryTrackerMap.end(); ++iter)
      {
        std::cout << iter->first << std::endl;
      }
    }
  }
  return;
}

std::vector<int> Fun4AllMemoryTracker::GetMemoryVector(const std::string &name) const
{
  std::vector<int> memvec;
  auto iter = mMemoryTrackerMap.find(name);
  if (iter != mMemoryTrackerMap.end())
  {
    memvec = iter->second;
  }
  return memvec;
}
