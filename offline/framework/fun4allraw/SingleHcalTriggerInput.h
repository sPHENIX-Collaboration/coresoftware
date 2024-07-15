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

 private:
  Packet **plist{nullptr};
};

#endif
