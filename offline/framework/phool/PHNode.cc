//  Author: Matthias Messer

#include "PHNode.h"

#include "phool.h"

#include <TSystem.h>

// boost stacktrace header causes a shadow warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

#include <iostream>

PHNode::PHNode(const std::string& n)
  : PHNode(n, "")
{
}

PHNode::PHNode(const std::string& n, const std::string& typ)
  :  objecttype(typ)
{
  int badnode = 0;
  if (n.find('.') != std::string::npos)
  {
    std::cout << PHWHERE << " No nodenames containing decimal point possible: "
         << n << std::endl;
    badnode = 1;
  }
  if (n.empty())
  {
    std::cout << PHWHERE << "Empty string as nodename given" << std::endl;
    badnode = 1;
  }
  if (n.find(' ') != std::string::npos)
  {
    badnode = 1;
    std::cout << PHWHERE << "No nodenames with spaces" << std::endl;
  }
  if (badnode)
  {
    std::cout << "Here is the stacktrace: " << std::endl;
    std::cout << boost::stacktrace::stacktrace();
    std::cout << "Check the stacktrace for the guilty party (typically #2)" << std::endl;
    gSystem->Exit(1);
  }
  name = n;
  return;
}

PHNode::~PHNode()
{
  if (parent)
  {
    parent->forgetMe(this);
  }
}

// Implementation of external functions.
std::ostream&
operator<<(std::ostream& stream, const PHNode& node)
{
  stream << node.getType() << " : " << node.getName() << " class " << node.getClass();

  return stream;
}
