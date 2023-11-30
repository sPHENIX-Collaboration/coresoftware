#ifndef FUN4ALLRAW_SINGLEZDCINPUT_H
#define FUN4ALLRAW_SINGLEZDCINPUT_H

#include "SinglePrdfInput.h"

#include <map>
#include <set>
#include <string>
#include <utility>  // for pair
#include <vector>

class Eventiterator;
class Fun4AllPrdfInputPoolManager;
class Packet;

class SingleZdcInput : public SinglePrdfInput
{
 public:
  explicit SingleZdcInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman);
  ~SingleZdcInput() override;
  void FillPool(const unsigned int nevents) override;

 private:
  int majority_eventnumber();
  int majority_beamclock();
  void adjust_eventnumber_offset(const int decided_evtno);
  struct PacketInfo
  {
    std::vector<Packet *> PacketVector;
    unsigned int EventFoundCounter = 0;
  };
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  int *m_PacketEventNumberOffset = nullptr;  // packet event counters start at 0 but we start with event number 1
  std::map<int, std::vector<Packet *>> m_PacketMap;
  std::set<int> m_EvtSet;
  std::vector<std::pair<int, int>> m_Event;
};

#endif
