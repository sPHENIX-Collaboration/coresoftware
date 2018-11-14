#pragma once

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpCaloTriggerInfo : public DumpObject
{
 public:
  DumpCaloTriggerInfo(const std::string &NodeName);
  virtual ~DumpCaloTriggerInfo() {}

 protected:
   int process_Node(PHNode *mynode);
};
