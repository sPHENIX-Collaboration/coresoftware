#ifndef NODEDUMP_DUMPGL1PACKET_H
#define NODEDUMP_DUMPGL1PACKET_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpGl1Packet : public DumpObject
{
 public:
  explicit DumpGl1Packet(const std::string &NodeName);
  ~DumpGl1Packet() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
