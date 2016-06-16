#ifndef DUMPPdbParameterMap_H__
#define DUMPPdbParameterMap_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPdbParameterMap : public DumpObject
{
 public:
  DumpPdbParameterMap(const std::string &NodeName);
  virtual ~DumpPdbParameterMap() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPPdbParameterMap_H__ */

