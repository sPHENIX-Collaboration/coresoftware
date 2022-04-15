#ifndef NODEDUMP_DUMPSVTXPHG4PARTICLEMAP_H
#define NODEDUMP_DUMPSVTXPHG4PARTICLEMAP_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpSvtxPHG4ParticleMap : public DumpObject
{
 public:
  DumpSvtxPHG4ParticleMap(const std::string &NodeName);
  ~DumpSvtxPHG4ParticleMap() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
