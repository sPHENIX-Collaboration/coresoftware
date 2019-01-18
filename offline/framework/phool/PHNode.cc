//  Author: Matthias Messer

#include "PHNode.h"

#include "phool.h"

#include <cstdlib>
#include <iostream>

using namespace std;

PHNode::PHNode()
  : parent(nullptr)
  , persistent(true)
  , type("PHNode")
  , reset_able(true)
{
  return;
}

PHNode::PHNode(const string& n)
  : parent(nullptr)
  , persistent(true)
  , type("PHNode")
  , reset_able(true)
{
  if (n.find(".") != string::npos)
  {
    cout << PHWHERE << " No nodenames containing decimal point possible: "
         << n << endl;
    exit(1);
  }
  name = n;
  return;
}

PHNode::PHNode(const string& n, const string& typ)
  : parent(nullptr)
  , persistent(true)
  , type("PHNode")
  , objecttype(typ)
  , reset_able(true)
{
  if (n.find(".") != string::npos)
  {
    cout << PHWHERE << " No nodenames containing decimal point possible: "
         << n << endl;
    exit(1);
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
