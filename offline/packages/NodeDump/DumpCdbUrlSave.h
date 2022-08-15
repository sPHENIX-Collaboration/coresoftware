#ifndef NODEDUMP_DUMPCDBURLSAVE_H
#define NODEDUMP_DUMPCDBURLSAVE_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpCdbUrlSave : public DumpObject
{
 public:
  explicit DumpCdbUrlSave(const std::string &NodeName);
  ~DumpCdbUrlSave() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
