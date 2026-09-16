#ifndef NODEDUMP_DUMPTPC_POLYTRACKCONTAINER_H
#define NODEDUMP_DUMPTPC_POLYTRACKCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpc_PolyTrackContainer : public DumpObject
{
 public:
  explicit DumpTpc_PolyTrackContainer(const std::string& NodeName);
  ~DumpTpc_PolyTrackContainer() override = default;

 protected:
  int process_Node(PHNode* myNode) override;
};

#endif
