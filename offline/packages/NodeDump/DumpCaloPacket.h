#ifndef NODEDUMP_DUMPCALOPACKET_H
#define NODEDUMP_DUMPCALOPACKET_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpCaloPacket : public DumpObject
{
 public:
  explicit DumpCaloPacket(const std::string &NodeName);
  ~DumpCaloPacket() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
