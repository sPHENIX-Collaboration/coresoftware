#ifndef NODEDUMP_DUMPPDBPARAMETERMAPCONTAINER_H
#define NODEDUMP_DUMPPDBPARAMETERMAPCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPdbParameterMapContainer : public DumpObject
{
 public:
  explicit DumpPdbParameterMapContainer(const std::string &NodeName);
  ~DumpPdbParameterMapContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
