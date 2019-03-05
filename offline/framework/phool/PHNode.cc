//  Author: Matthias Messer

#include "PHNode.h"

#include "phool.h"

#include <TSystem.h>

#include <boost/stacktrace.hpp>

#include <cstdlib>
#include <iostream>

using namespace std;

PHNode::PHNode(const string& n): PHNode(n,"")
{}

PHNode::PHNode(const string& n, const string& typ)
  : parent(nullptr)
  , persistent(true)
  , type("PHNode")
  , objecttype(typ)
  , reset_able(true)
{
  int badnode = 0;
  if (n.find(".") != string::npos)
  {
    cout << PHWHERE << " No nodenames containing decimal point possible: "
         << n << endl;
    badnode = 1;
  }
  if (n.empty())
  {
    cout << PHWHERE << "Empty string as nodename given" << endl;
    badnode = 1;
  }
  if (n.find(" ") != string::npos)
  {
    badnode = 1;
    cout << PHWHERE << "No nodenames with spaces" << endl;
  }
  if (badnode)
  {
    cout << "Here is the stacktrace: " << endl;
    cout << boost::stacktrace::stacktrace();
    cout << "Check the stacktrace for the guilty party (typically #2)" << endl;
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
