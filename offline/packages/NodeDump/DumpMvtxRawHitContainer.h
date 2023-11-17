#ifndef NODEDUMP_DUMPMVTXRAWHITCONTAINER_H
#define NODEDUMP_DUMPMVTXRAWHITCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMvtxRawHitContainer : public DumpObject
{
 public:
  explicit DumpMvtxRawHitContainer(const std::string &NodeName);
  ~DumpMvtxRawHitContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
