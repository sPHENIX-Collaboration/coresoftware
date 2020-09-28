#ifndef NODEDUMP_DUMPPHHEPMCGENEVENTMAP_H
#define NODEDUMP_DUMPPHHEPMCGENEVENTMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHHepMCGenEventMap : public DumpObject
{
 public:
  DumpPHHepMCGenEventMap(const std::string &NodeName);
  virtual ~DumpPHHepMCGenEventMap() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
