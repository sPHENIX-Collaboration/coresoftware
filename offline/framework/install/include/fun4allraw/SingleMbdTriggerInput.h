#ifndef FUN4ALLRAW_SINGLEMBDTRIGGERINPUT_H
#define FUN4ALLRAW_SINGLEMBDTRIGGERINPUT_H

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

class SingleMbdTriggerInput : public SingleTriggerInput
{
 public:
  explicit SingleMbdTriggerInput(const std::string &name);
  ~SingleMbdTriggerInput() override;
  void FillPool(const unsigned int) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

 private:
  Packet **plist{nullptr};
};

#endif
