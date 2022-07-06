#ifndef NODEDUMP_DUMPPARTICLEFLOWELEMENTCONTAINER_H
#define NODEDUMP_DUMPPARTICLEFLOWELEMENTCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpParticleFlowElementContainer : public DumpObject
{
 public:
  explicit DumpParticleFlowElementContainer(const std::string &NodeName);
  ~DumpParticleFlowElementContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
