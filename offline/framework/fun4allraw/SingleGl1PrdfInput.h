#ifndef FUN4ALLRAW_SINGLEGL1PRDFINPUT_H
#define FUN4ALLRAW_SINGLEGL1PRDFINPUT_H

#include "SinglePrdfInput.h"

#include <array>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>  // for pair
#include <vector>

class Eventiterator;
class Fun4AllPrdfInputPoolManager;
class Fun4AllPrdfInputTriggerManager;
class Packet;

class SingleGl1PrdfInput : public SinglePrdfInput
{
 public:
  explicit SingleGl1PrdfInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman);
  explicit SingleGl1PrdfInput(const std::string &name, Fun4AllPrdfInputTriggerManager *inman);
  ~SingleGl1PrdfInput() override;
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
  std::map<uint64_t, std::vector<Packet *>> m_PacketMap;
  std::set<int> m_EvtSet;
  std::vector<std::pair<int, uint64_t>> m_Event;
};

#endif
