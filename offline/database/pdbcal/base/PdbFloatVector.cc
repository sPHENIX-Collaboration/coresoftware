#include "PdbFloatVector.hh"
#include <iostream>

using namespace std;

PdbFloatVector::~PdbFloatVector()
{
  _floatvec.clear();
  return;
}

void PdbFloatVector::print() const
{
  for( vector<float>::const_iterator iter = _floatvec.begin(); iter != _floatvec.end(); iter++ )
    { cout << *iter << endl; }
}

void PdbFloatVector::add_float(float fl){
  _floatvec.push_back(fl);
}

float PdbFloatVector::getValue(int pos) const {
  return _floatvec.at(pos);
}

