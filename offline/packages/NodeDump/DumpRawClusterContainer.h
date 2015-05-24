#ifndef DUMPRAWCLUSTERCONTAINER_H__
#define DUMPRAWCLUSTERCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawClusterContainer : public DumpObject
{
 public:
  DumpRawClusterContainer(const std::string &NodeName);
  virtual ~DumpRawClusterContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPRAWCLUSTERCONTAINER_H__ */

