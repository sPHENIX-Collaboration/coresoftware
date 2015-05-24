#ifndef DUMPPHG4INEVENT_H__
#define DUMPPHG4INEVENT_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPHG4InEvent : public DumpObject
{
 public:
  DumpPHG4InEvent(const std::string &NodeName);
  virtual ~DumpPHG4InEvent() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* __DUMPPHG4INEVENT_H__ */

