#ifndef NODEDUMP_DUMPPHG4SCINTILLATORSLATCONTAINER_H
#define NODEDUMP_DUMPPHG4SCINTILLATORSLATCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4ScintillatorSlatContainer : public DumpObject
{
 public:
  explicit DumpPHG4ScintillatorSlatContainer(const std::string &NodeName);
  ~DumpPHG4ScintillatorSlatContainer() override {}

 protected:
  int process_Node(PHNode *mynode) override;
};

#endif
