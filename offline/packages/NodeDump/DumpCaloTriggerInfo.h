#pragma once

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpCaloTriggerInfo : public DumpObject
{
 public:
  explicit DumpCaloTriggerInfo(const std::string &NodeName);
  ~DumpCaloTriggerInfo() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};
