#ifndef NODEDUMP_DUMPINTTDEADMAP_H
#define NODEDUMP_DUMPINTTDEADMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpInttDeadMap : public DumpObject
{
 public:
  explicit DumpInttDeadMap(const std::string &NodeName);
  ~DumpInttDeadMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
