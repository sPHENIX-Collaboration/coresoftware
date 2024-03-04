#ifndef NODEDUMP_DUMPJETCONTAINER_H
#define NODEDUMP_DUMPJETCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpJetContainer : public DumpObject
{
 public:
  explicit DumpJetContainer(const std::string &NodeName);
  ~DumpJetContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
