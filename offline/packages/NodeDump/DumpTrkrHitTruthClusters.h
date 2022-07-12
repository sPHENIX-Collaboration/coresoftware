#ifndef NODEDUMP_DUMPTRKRHITTRUTHCLUSTERS_H
#define NODEDUMP_DUMPTRKRHITTRUTHCLUSTERS_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrHitTruthClusters : public DumpObject
{
 public:
  explicit DumpTrkrHitTruthClusters(const std::string &NodeName);
  ~DumpTrkrHitTruthClusters() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
