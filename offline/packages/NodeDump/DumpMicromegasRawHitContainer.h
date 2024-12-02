#ifndef NODEDUMP_DUMPMICROMEGASRAWHITCONTAINER_H
#define NODEDUMP_DUMPMICROMEGASRAWHITCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMicromegasRawHitContainer : public DumpObject
{
 public:
  explicit DumpMicromegasRawHitContainer(const std::string &NodeName);
  ~DumpMicromegasRawHitContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
