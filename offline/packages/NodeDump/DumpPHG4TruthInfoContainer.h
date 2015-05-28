#ifndef DUMPPHG4TRUTHINFOCONTAINER_H__
#define DUMPPHG4TRUTHINFOCONTAINER_H__

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

#endif /* __DUMPPHG4TRUTHCONTAINER_H__ */

