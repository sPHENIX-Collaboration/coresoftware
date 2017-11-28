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

PHNode::PHNode(const PHNode& phn)
  : parent(nullptr)
  , persistent(phn.persistent)
  , type(phn.type)
  , objecttype(phn.objecttype)
  , name(phn.name)
  , objectclass(phn.objectclass)
  , reset_able(phn.reset_able)
{
  cout << "copy ctor not implemented because of pointer to parent" << endl;
  cout << "which needs implementing for this to be reasonable" << endl;
  exit(1);
}

PHNode&
PHNode::operator=(const PHNode&)
{
  cout << "= operator not implemented because of pointer to parent" << endl;
  cout << "which needs implementing for this to be reasonable" << endl;
  exit(1);
}

// void
// PHNode::setResetFlag(const int val)
// {
//   reset_able = (val) ? true : false;
// }

// Implementation of external functions.
std::ostream&
operator<<(std::ostream& stream, const PHNode& node)
{
  stream << node.getType() << " : " << node.getName() << " class " << node.getClass();

  return stream;
}
