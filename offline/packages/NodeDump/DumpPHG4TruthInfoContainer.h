#ifndef NODEDUMP_DUMPPHG4TRUTHINFOCONTAINER_H
#define NODEDUMP_DUMPPHG4TRUTHINFOCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4TruthInfoContainer : public DumpObject
{
 public:
  DumpPHG4TruthInfoContainer(const std::string &NodeName);
  virtual ~DumpPHG4TruthInfoContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
