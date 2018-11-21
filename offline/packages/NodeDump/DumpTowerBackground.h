#pragma once

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpTowerBackground : public DumpObject
{
 public:
  DumpTowerBackground(const std::string &NodeName);
  virtual ~DumpTowerBackground() {}

 protected:
   int process_Node(PHNode *mynode);
};

