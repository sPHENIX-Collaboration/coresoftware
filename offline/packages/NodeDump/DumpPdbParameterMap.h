#ifndef NODEDUMP_DUMPPDBPARAMETERMAP_H
#define NODEDUMP_DUMPPDBPARAMETERMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPdbParameterMap : public DumpObject
{
 public:
  DumpPdbParameterMap(const std::string &NodeName);
  virtual ~DumpPdbParameterMap() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
