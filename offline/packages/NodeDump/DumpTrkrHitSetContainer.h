#ifndef NODEDUMP_DUMPTRKRHITSETCONTAINER_H
#define NODEDUMP_DUMPTRKRHITSETCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrHitSetContainer : public DumpObject
{
 public:
  DumpTrkrHitSetContainer(const std::string &NodeName);
  virtual ~DumpTrkrHitSetContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
