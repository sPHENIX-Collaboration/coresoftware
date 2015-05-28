#ifndef DUMPRAWTOWERCONTAINER_H__
#define DUMPRAWTOWERCONTAINER_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpRawTowerContainer : public DumpObject
{
 public:
  DumpRawTowerContainer(const std::string &NodeName);
  virtual ~DumpRawTowerContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPRAWTOWERCONTAINER_H__ */

