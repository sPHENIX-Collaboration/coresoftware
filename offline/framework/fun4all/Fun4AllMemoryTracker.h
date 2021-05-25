// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLMEMORYTRACKER_H
#define FUN4ALL_FUN4ALLMEMORYTRACKER_H

#include "Fun4AllBase.h"

#include <map>
#include <string>
#include <vector>

class Fun4AllMemoryTracker : public Fun4AllBase
{
 public:
  static Fun4AllMemoryTracker *instance()
  {
    if (mInstance) return mInstance;
    mInstance = new Fun4AllMemoryTracker();
    return mInstance;
  }
  ~Fun4AllMemoryTracker() override;
  void Snapshot(const std::string &trackername, const std::string &group = "");
  void Start(const std::string &trackername, const std::string &group = "");
  void Stop(const std::string &trackername, const std::string &group = "");

  int GetRSSMemory() const;
  void PrintMemoryTracker(const std::string &name = "") const;
  std::vector<int> GetMemoryVector(const std::string &name) const;

 private:
  Fun4AllMemoryTracker();
  std::string CreateFullTrackerName(const std::string &trackername, const std::string &group = "");
  static Fun4AllMemoryTracker *mInstance;
  std::map<std::string, std::vector<int>> mMemoryTrackerMap;
  std::map<std::string, int> mStartMem;
};

#endif
