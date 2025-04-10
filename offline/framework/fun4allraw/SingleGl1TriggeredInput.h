#ifndef FUN4ALLRAW_SINGLEGL1TRIGGEREDINPUT_H
#define FUN4ALLRAW_SINGLEGL1TRIGGEREDINPUT_H

#include "SingleTriggeredInput.h"

#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

class OfflinePacket;
class Packet;
class PHCompositeNode;

class SingleGl1TriggeredInput : public SingleTriggeredInput
{
 public:
  explicit SingleGl1TriggeredInput(const std::string &name);
  ~SingleGl1TriggeredInput() override = default;
  ;
  void FillPool(const unsigned int) override;
  // void CleanupUsedPackets(const int eventno);
  // void ClearCurrentEvent();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  uint64_t GetClock(Event *evt) override;

 private:
};

#endif
