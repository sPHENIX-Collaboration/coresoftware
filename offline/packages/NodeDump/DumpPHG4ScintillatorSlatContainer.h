#ifndef DumpPHG4ScintillatorSlatContainer_H__
#define DumpPHG4ScintillatorSlatContainer_H__

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

#endif /* DUMPRAWCLUSTERCONTAINER_H__ */

