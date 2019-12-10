#ifndef NODEDUMP_DUMPPDBPARAMETERMAPCONTAINER_H
#define NODEDUMP_DUMPPDBPARAMETERMAPCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPdbParameterMapContainer : public DumpObject
{
 public:
  DumpPdbParameterMapContainer(const std::string &NodeName);
  virtual ~DumpPdbParameterMapContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
