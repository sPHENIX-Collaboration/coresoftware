// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLMONITORING_H
#define FUN4ALL_FUN4ALLMONITORING_H

#include "Fun4AllBase.h"

#include <cstdint>
#include <string>

class Fun4AllMonitoring : public Fun4AllBase
{
 public:
  static Fun4AllMonitoring *instance()
  {
    if (mInstance) return mInstance;
    mInstance = new Fun4AllMonitoring();
    return mInstance;
  }
  ~Fun4AllMonitoring() override = default;
  void Snapshot(const std::string &what = "AfterProcessEvent");

  void PrintsMaps() const;

  void Get_Memory();
  void OutFileName(const std::string &fname);

 private:
  Fun4AllMonitoring();
  static Fun4AllMonitoring *mInstance;
  uint64_t mEvent = 0;
  uint64_t mHeapPss = 0;
  uint64_t mMMapPSS = 0;
  uint64_t mOtherPss = 0;
  std::string mOutFileName;
};

#endif
