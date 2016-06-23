#ifndef DUMPPHG4BLOCKGEOMCONTAINER_H__
#define DUMPPHG4BLOCKGEOMCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4BlockGeomContainer : public DumpObject
{
 public:
  DumpPHG4BlockGeomContainer(const std::string &NodeName);
  virtual ~DumpPHG4BlockGeomContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPPHG4BLOCKGEOMCONTAINER_H__ */

