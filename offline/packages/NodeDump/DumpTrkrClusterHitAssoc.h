#ifndef NODEDUMP_DUMPTRKRCLUSTERHITASSOC_H
#define NODEDUMP_DUMPTRKRCLUSTERHITASSOC_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTrkrClusterHitAssoc : public DumpObject
{
 public:
  DumpTrkrClusterHitAssoc(const std::string &NodeName);
  virtual ~DumpTrkrClusterHitAssoc() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
