#ifndef NODEDUMP_DUMPTPCCROSSINGDECISIONCONTAINER_H
#define NODEDUMP_DUMPTPCCROSSINGDECISIONCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTpcCrossingDecisionContainer : public DumpObject
{
 public:
  explicit DumpTpcCrossingDecisionContainer(const std::string& NodeName);
  ~DumpTpcCrossingDecisionContainer() override = default;

 protected:
  int process_Node(PHNode* myNode) override;
};

#endif
