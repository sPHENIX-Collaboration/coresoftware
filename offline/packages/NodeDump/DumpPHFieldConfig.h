#ifndef NODEDUMP_DUMPPHFIELDCONFIG_H
#define NODEDUMP_DUMPPHFIELDCONFIG_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHFieldConfig : public DumpObject
{
 public:
  explicit DumpPHFieldConfig(const std::string &NodeName);
  ~DumpPHFieldConfig() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
