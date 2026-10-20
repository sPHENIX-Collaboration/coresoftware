#ifndef NODEDUMP_DUMPLASERCLUSTERCONTAINER_H
#define NODEDUMP_DUMPLASERCLUSTERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpLaserClusterContainer : public DumpObject
{
 public:
  explicit DumpLaserClusterContainer(const std::string& NodeName);
  ~DumpLaserClusterContainer() override = default;

 protected:
  int process_Node(PHNode* mynode) override;
};

#endif
