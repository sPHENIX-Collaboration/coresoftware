#ifndef NODEDUMP_DUMPPDBPARAMETERMAP_H
#define NODEDUMP_DUMPPDBPARAMETERMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPdbParameterMap : public DumpObject
{
 public:
  explicit DumpPdbParameterMap(const std::string &NodeName);
  ~DumpPdbParameterMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
