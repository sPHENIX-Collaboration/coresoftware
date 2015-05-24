#ifndef __DUMPRUNHEADER_H__
#define __DUMPRUNHEADER_H__

#include <DumpObject.h>

#include <string>

class PHNode;

class DumpRunHeader : public DumpObject
{
 public:
  DumpRunHeader(const std::string &NodeName);
  virtual ~DumpRunHeader() {}

 protected:
   int process_Node(PHNode *mynode);
   int node_written;
};

#endif /* __DUMPRUNHEADER_H__ */

