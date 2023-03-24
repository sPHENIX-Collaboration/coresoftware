#ifndef NODEDUMP_DUMPTOWERINFOCONTAINER_H
#define NODEDUMP_DUMPTOWERINFOCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTowerInfoContainer : public DumpObject
{
 public:
  explicit DumpTowerInfoContainer(const std::string &NodeName);
  ~DumpTowerInfoContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
