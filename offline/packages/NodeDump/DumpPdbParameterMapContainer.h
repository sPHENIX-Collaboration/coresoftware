#ifndef DUMPPdbParameterMapContainer_H__
#define DUMPPdbParameterMapContainer_H__

#include "DumpObject.h"

#include <string>

class PHNode;

class DumpPdbParameterMapContainer : public DumpObject
{
 public:
  DumpPdbParameterMapContainer(const std::string &NodeName);
  virtual ~DumpPdbParameterMapContainer() {}

 protected:
   int process_Node(PHNode *mynode);
};

#endif /* DUMPPdbParameterMapContainer_H__ */

