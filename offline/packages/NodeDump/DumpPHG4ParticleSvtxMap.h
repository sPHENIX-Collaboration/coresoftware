#ifndef NODEDUMP_PHG4PARTICLESVTXMAP_H
#define NODEDUMP_PHG4PARTICLESVTXMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4ParticleSvtxMap : public DumpObject
{
 public:
  explicit DumpPHG4ParticleSvtxMap(const std::string &NodeName);
  ~DumpPHG4ParticleSvtxMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
