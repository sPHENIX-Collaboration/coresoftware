#pragma once

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTowerBackground : public DumpObject
{
 public:
  explicit DumpTowerBackground(const std::string &NodeName);
  ~DumpTowerBackground() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};
