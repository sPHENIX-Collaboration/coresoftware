#ifndef NODEDUMP_DUMPTPC_POLYCLUSTERCONTAINER_H
#define NODEDUMP_DUMPTPC_POLYCLUSTERCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpc_PolyClusterContainer : public DumpObject
{
 public:
  explicit DumpTpc_PolyClusterContainer(const std::string& NodeName);
  ~DumpTpc_PolyClusterContainer() override = default;

 protected:
  int process_Node(PHNode* myNode) override;
};

#endif
