#ifndef NODEDUMP_DUMPPHG4INEVENT_H
#define NODEDUMP_DUMPPHG4INEVENT_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4InEvent : public DumpObject
{
 public:
  explicit DumpPHG4InEvent(const std::string &NodeName);
  ~DumpPHG4InEvent() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
