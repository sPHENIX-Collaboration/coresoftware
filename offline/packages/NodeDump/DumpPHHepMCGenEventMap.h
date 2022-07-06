#ifndef NODEDUMP_DUMPPHHEPMCGENEVENTMAP_H
#define NODEDUMP_DUMPPHHEPMCGENEVENTMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHHepMCGenEventMap : public DumpObject
{
 public:
  explicit DumpPHHepMCGenEventMap(const std::string &NodeName);
  ~DumpPHHepMCGenEventMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
