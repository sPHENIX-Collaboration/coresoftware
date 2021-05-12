#ifndef NODEDUMP_DUMPTRKRHITTRUTHASSOC_H
#define NODEDUMP_DUMPTRKRHITTRUTHASSOC_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrHitTruthAssoc : public DumpObject
{
 public:
  DumpTrkrHitTruthAssoc(const std::string &NodeName);
  virtual ~DumpTrkrHitTruthAssoc() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
