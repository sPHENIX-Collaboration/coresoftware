#ifndef NODEDUMP_DUMPTRKRCLUSTERCROSSINGASSOC_H
#define NODEDUMP_DUMPTRKRCLUSTERCROSSINGASSOC_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrClusterCrossingAssoc : public DumpObject
{
 public:
  explicit DumpTrkrClusterCrossingAssoc(const std::string &NodeName);
  ~DumpTrkrClusterCrossingAssoc() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
