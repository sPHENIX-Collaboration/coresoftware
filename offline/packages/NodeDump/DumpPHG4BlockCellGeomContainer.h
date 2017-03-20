#ifndef DUMPPHG4BLOCKCELLGEOMCONTAINER_H__
#define DUMPPHG4BLOCKCELLGEOMCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4BlockCellGeomContainer : public DumpObject
{
 public:
  DumpPHG4BlockCellGeomContainer(const std::string &NodeName);
  virtual ~DumpPHG4BlockCellGeomContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPPHG4BLOCKCELLGEOMCONTAINER_H__ */

