#ifndef FUN4ALLRAW_SINGLEHCALTRIGGERINPUT_H
#define FUN4ALLRAW_SINGLEHCALTRIGGERINPUT_H

#include "SingleTriggerInput.h"

#include <cstdint>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class OfflinePacket;
class Packet;
class PHCompositeNode;

class SingleHcalTriggerInput : public SingleTriggerInput
{
 public:
  explicit SingleHcalTriggerInput(const std::string &name);
  ~SingleHcalTriggerInput() override;
  void FillPool(const unsigned int keep) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void CleanupUsedLocalPackets(const int eventno);
  void SetFEMClockProblemFlag(bool b = true) { m_FEMClockProblemFlag = b; }
  bool FEMClockProblemFlag() const { return m_FEMClockProblemFlag; }
  void SetClockReferencePacket(const int i) { m_ClockReferencePacket = i; }
  int ClockReferencePacket() const { return m_ClockReferencePacket; }

 private:
  void CheckFEMClock();
  void CheckFEMEventNumber();
  int ShiftEvents(int pktid, int offset);
  int m_ClockReferencePacket{0};
  bool m_FEMClockProblemFlag{false};
  Packet **plist{nullptr};
  std::set<int> m_BadBCOPacketSet;
  std::map<int, std::vector<OfflinePacket *>> m_LocalPacketMap;
  std::map<int, uint64_t> m_EventRefBCO;
};

#endif
