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
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

 private:
  Packet **plist{nullptr};
};

#endif
