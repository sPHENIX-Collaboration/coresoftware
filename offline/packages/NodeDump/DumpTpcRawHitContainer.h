#ifndef NODEDUMP_DUMPTPCRAWHITCONTAINER_H
#define NODEDUMP_DUMPTPCRAWHITCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpcRawHitContainer : public DumpObject
{
 public:
  explicit DumpTpcRawHitContainer(const std::string &NodeName);
  ~DumpTpcRawHitContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
