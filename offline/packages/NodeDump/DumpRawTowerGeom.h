#ifndef DUMPRAWTOWERGEOM_H__
#define DUMPRAWTOWERGEOM_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerGeom : public DumpObject
{
 public:
  DumpRawTowerGeom(const std::string &NodeName);
  virtual ~DumpRawTowerGeom() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPRAWTOWERGEOM_H__ */

