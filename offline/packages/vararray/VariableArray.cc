#include "VariableArray.h"

#include <ostream>

using namespace std;

VariableArray::VariableArray(const unsigned int idval)
  : id(idval)
  , nVal(0)
  , sval(nullptr)
{
  return;
}

VariableArray::~VariableArray()
{
  Reset();
  return;
}

void VariableArray::identify(ostream &os) const
{
  os << "contain " << nVal << " values" << endl;
  for (unsigned int i = 0; i < nVal; i++)
  {
    os << "n: " << i << " val: " << sval[i] << endl;
  }
  return;
}

void VariableArray::set_val(const vector<short> &vec)
{
  nVal = vec.size();
  sval = new short[nVal];
  vector<short>::const_iterator iter;
  unsigned int i = 0;
  for (iter = vec.begin(); iter != vec.end(); ++iter)
  {
    sval[i++] = *iter;
  }
  return;
}

void VariableArray::Reset()
{
  delete[] sval;
  sval = nullptr;
  nVal = 0;
  return;
}
