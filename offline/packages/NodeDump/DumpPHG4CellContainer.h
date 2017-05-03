#ifndef DUMPPHG4CELLCONTAINER_H__
#define DUMPPHG4CELLCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4CellContainer : public DumpObject
{
 public:
  DumpPHG4CellContainer(const std::string &NodeName);
  virtual ~DumpPHG4CellContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* __DUMPPHG4CELLCONTAINER_H__ */

