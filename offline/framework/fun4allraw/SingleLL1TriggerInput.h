#ifndef FUN4ALLRAW_SINGLELL1TRIGGERINPUT_H
#define FUN4ALLRAW_SINGLELL1TRIGGERINPUT_H

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

class SingleLL1TriggerInput : public SingleTriggerInput
{
 public:
  explicit SingleLL1TriggerInput(const std::string &name);
  ~SingleLL1TriggerInput() override;
  void FillPool(const unsigned int keep) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents(const unsigned int keep);
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  //  void ConfigureStreamingInputManager() override;

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};

  std::set<int> m_EventNumber;
  std::map<int, std::vector<OfflinePacket *>> m_LL1PacketMap;
  std::set<int> m_EventStack;
};

#endif
