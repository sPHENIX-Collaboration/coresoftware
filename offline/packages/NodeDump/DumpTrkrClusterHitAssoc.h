#ifndef NODEDUMP_DUMPTRKRCLUSTERHITASSOC_H
#define NODEDUMP_DUMPTRKRCLUSTERHITASSOC_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrClusterHitAssoc : public DumpObject
{
 public:
  explicit DumpTrkrClusterHitAssoc(const std::string &NodeName);
  ~DumpTrkrClusterHitAssoc() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
