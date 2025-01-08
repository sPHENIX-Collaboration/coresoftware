#ifndef NODEDUMP_DUMPMVTXRAWEVTHEADER_H
#define NODEDUMP_DUMPMVTXRAWEVTHEADER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpMvtxRawEvtHeader : public DumpObject
{
 public:
  explicit DumpMvtxRawEvtHeader(const std::string &NodeName);
  ~DumpMvtxRawEvtHeader() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
