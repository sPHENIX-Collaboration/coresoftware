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
  void FillPool(const unsigned int keep) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

 private:
  unsigned int m_Gl1PacketNumberEventNumberDiff{0};
};

#endif
