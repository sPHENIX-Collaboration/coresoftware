#ifndef NODEDUMP_DUMPPHG4SCINTILLATORSLATCONTAINER_H
#define NODEDUMP_DUMPPHG4SCINTILLATORSLATCONTAINER_H

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4ScintillatorSlatContainer : public DumpObject
{
 public:
  DumpPHG4ScintillatorSlatContainer(const std::string &NodeName);
  virtual ~DumpPHG4ScintillatorSlatContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif

