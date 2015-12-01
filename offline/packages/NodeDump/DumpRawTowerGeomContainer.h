#ifndef DumpRawTowerGeomContainer_H__
#define DumpRawTowerGeomContainer_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerGeomContainer : public DumpObject
{
 public:
  DumpRawTowerGeomContainer(const std::string &NodeName);
  virtual ~DumpRawTowerGeomContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DumpRawTowerGeomContainer_H__ */

