#ifndef NODEDUMP_DUMPTRKRHITSETCONTAINER_H
#define NODEDUMP_DUMPTRKRHITSETCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrHitSetContainer : public DumpObject
{
 public:
  explicit DumpTrkrHitSetContainer(const std::string &NodeName);
  ~DumpTrkrHitSetContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
