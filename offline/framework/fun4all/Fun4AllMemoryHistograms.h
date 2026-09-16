// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLMEMORYHISTOGRAMS_H
#define FUN4ALL_FUN4ALLMEMORYHISTOGRAMS_H

#include "Fun4AllBase.h"


#include <map>
#include <string>
#include <vector>

class TH1;

class Fun4AllMemoryHistograms : public Fun4AllBase
{
 public:
  static Fun4AllMemoryHistograms *instance()
  {
    if (mInstance) return mInstance;
    mInstance = new Fun4AllMemoryHistograms();
    return mInstance;
  }
  ~Fun4AllMemoryHistograms() override = default;
  void Enable() {m_Enabled = true;}
  void CreateHistos();
  void FillHistos();
  void SaveHistos(const std::string &fname);
  void increase_TH1_channelcount();
  
 private:
  Fun4AllMemoryHistograms();
  static Fun4AllMemoryHistograms *mInstance;
  
  int m_EventCount {0};
  int m_NumChannels {10000};
  bool m_Enabled {false};
  bool m_HistosCreated {false};

  std::vector<TH1 *> m_histovector;
};

#endif
