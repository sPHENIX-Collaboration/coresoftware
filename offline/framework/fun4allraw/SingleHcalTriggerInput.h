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
  void FillPool(const unsigned int) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  //  void ConfigureStreamingInputManager() override;

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};

  std::set<int> m_EventNumber;
  std::map<int, std::vector<OfflinePacket *>> m_HcalPacketMap;
  std::set<int> m_EventStack;
};

#endif
