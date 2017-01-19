#ifndef DUMPJETMAP_H__
#define DUMPJETMAP_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpJetMap : public DumpObject
{
 public:
  DumpJetMap(const std::string &NodeName);
  virtual ~DumpJetMap() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* __DUMPJETMAP_H__ */

