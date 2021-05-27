#ifndef NODEDUMP_DUMPJETMAP_H
#define NODEDUMP_DUMPJETMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpJetMap : public DumpObject
{
 public:
  DumpJetMap(const std::string &NodeName);
  ~DumpJetMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
