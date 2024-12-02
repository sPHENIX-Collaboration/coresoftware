#ifndef NODEDUMP_DUMPZDCINFO_H
#define NODEDUMP_DUMPZDCINFO_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpZdcinfo : public DumpObject
{
 public:
  explicit DumpZdcinfo(const std::string &NodeName);
  virtual ~DumpZdcinfo() {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
