#ifndef NODEDUMP_DUMPTRKRCLUSTERCONTAINER_H
#define NODEDUMP_DUMPTRKRCLUSTERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrClusterContainer : public DumpObject
{
 public:
  explicit DumpTrkrClusterContainer(const std::string &NodeName);
  ~DumpTrkrClusterContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
