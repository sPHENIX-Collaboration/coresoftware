#ifndef FUN4ALLRAW_SINGLECEMCTRIGGERINPUT_H
#define FUN4ALLRAW_SINGLECEMCTRIGGERINPUT_H

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

class SingleCemcTriggerInput : public SingleTriggerInput
{
 public:
  explicit SingleCemcTriggerInput(const std::string &name);
  ~SingleCemcTriggerInput() override;
  void FillPool(const unsigned int keep) override;
  void CleanupUsedPackets(const int eventno) override;
  void ClearCurrentEvent() override;
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  int CheckFEMClocks();
  void CleanupUsedLocalPackets(const int eventno);
 
 private:
  Packet **plist{nullptr};
  std::map<int, std::vector<OfflinePacket *>> m_LocalPacketMap;
};

#endif
