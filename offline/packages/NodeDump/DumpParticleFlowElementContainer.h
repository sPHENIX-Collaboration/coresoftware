#ifndef NODEDUMP_DUMPPARTICLEFLOWELEMENTCONTAINER_H
#define NODEDUMP_DUMPPARTICLEFLOWELEMENTCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpParticleFlowElementContainer : public DumpObject
{
 public:
  DumpParticleFlowElementContainer(const std::string &NodeName);
  virtual ~DumpParticleFlowElementContainer() {}

 protected:
  int process_Node(PHNode *mynode);
};

#endif
