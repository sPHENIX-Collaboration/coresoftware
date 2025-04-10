#ifndef FUN4ALLRAW_SINGLEZDCTRIGGEREDINPUT_H
#define FUN4ALLRAW_SINGLEZDCTRIGGEREDINPUT_H

#include "SingleTriggeredInput.h"

#include <array>
#include <cstdint>  // for uint64_t
#include <deque>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

class Event;
class Eventiterator;
class Fun4AllPrdfInputTriggerManager;
class OfflinePacket;
class Packet;
class PHCompositeNode;

class SingleZdcTriggeredInput : public SingleTriggeredInput
{
 public:
  explicit SingleZdcTriggeredInput(const std::string &name);
  ~SingleZdcTriggeredInput() override = default;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void AddPacket(PHCompositeNode *topNode, OfflinePacket *newhit) override;

 private:
};

#endif
