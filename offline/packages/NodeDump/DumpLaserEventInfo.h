#ifndef NODEDUMP_DUMPLASEREVENTINFO_H
#define NODEDUMP_DUMPLASEREVENTINFO_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpLaserEventInfo : public DumpObject
{
 public:
  explicit DumpLaserEventInfo(const std::string& NodeName);
  ~DumpLaserEventInfo() override = default;

 protected:
  int process_Node(PHNode* myNode) override;
};

#endif
