#ifndef NODEDUMP_DUMPINTTRAWHITCONTAINER_H
#define NODEDUMP_DUMPINTTRAWHITCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpInttRawHitContainer : public DumpObject
{
 public:
  explicit DumpInttRawHitContainer(const std::string &NodeName);
  ~DumpInttRawHitContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
