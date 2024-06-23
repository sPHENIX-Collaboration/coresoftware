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
  bool GetSomeMoreEvents(const unsigned int keep);
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

 private:

  std::set<int> m_EventNumber;
  std::set<int> m_EventStack;
};

#endif
