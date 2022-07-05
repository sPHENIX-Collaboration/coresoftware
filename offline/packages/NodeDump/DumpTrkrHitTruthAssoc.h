#ifndef NODEDUMP_DUMPTRKRHITTRUTHASSOC_H
#define NODEDUMP_DUMPTRKRHITTRUTHASSOC_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrHitTruthAssoc : public DumpObject
{
 public:
  explicit DumpTrkrHitTruthAssoc(const std::string &NodeName);
  ~DumpTrkrHitTruthAssoc() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
