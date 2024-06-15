#ifndef NODEDUMP_DUMPCALOPACKETCONTAINER_H
#define NODEDUMP_DUMPCALOPACKETCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpCaloPacketContainer : public DumpObject
{
 public:
  explicit DumpCaloPacketContainer(const std::string &NodeName);
  ~DumpCaloPacketContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
