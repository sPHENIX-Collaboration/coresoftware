#ifndef NODEDUMP_DUMPRUNHEADER_H
#define NODEDUMP_DUMPRUNHEADER_H

#include <DumpObject.h>

#include <string>

class PHNode;

class DumpRunHeader : public DumpObject
{
 public:
  DumpRunHeader(const std::string &NodeName);
  ~DumpRunHeader() override {}

 protected:
  int process_Node(PHNode *mynode) override;
  int node_written = 0;
};

#endif
