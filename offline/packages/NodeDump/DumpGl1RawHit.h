#ifndef NODEDUMP_DUMPGL1RAWHIT_H
#define NODEDUMP_DUMPGL1RAWHIT_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpGl1RawHit : public DumpObject
{
 public:
  explicit DumpGl1RawHit(const std::string &NodeName);
  ~DumpGl1RawHit() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
