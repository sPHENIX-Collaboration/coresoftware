#ifndef NODEDUMP_DUMPTPC_POLYTRACKVERTEXCONTAINER_H
#define NODEDUMP_DUMPTPC_POLYTRACKVERTEXCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpc_PolyTrackVertexContainer : public DumpObject
{
 public:
  explicit DumpTpc_PolyTrackVertexContainer(const std::string& NodeName);
  ~DumpTpc_PolyTrackVertexContainer() override = default;

 protected:
  int process_Node(PHNode* myNode) override;
};

#endif
