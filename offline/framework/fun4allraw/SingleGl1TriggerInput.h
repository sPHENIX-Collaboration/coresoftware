#ifndef FUN4ALLRAW_SINGLEGL1TRIGGERINPUT_H
#define FUN4ALLRAW_SINGLEGL1TRIGGERINPUT_H

#include "SingleTriggerInput.h"

#include <cstdint>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class OfflinePacket;
class PHCompositeNode;

class SingleGl1TriggerInput : public SingleTriggerInput
{
 public:
  explicit SingleGl1TriggerInput(const std::string &name);
  ~SingleGl1TriggerInput() override;
  void FillPool(const unsigned int) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  //  void ConfigureStreamingInputManager() override;

 private:
  unsigned int m_NumSpecialEvents{0};

  std::set<int> m_EventNumber;
  std::map<int, std::vector<OfflinePacket *>> m_Gl1PacketMap;
  std::set<int> m_EventStack;
};

#endif
