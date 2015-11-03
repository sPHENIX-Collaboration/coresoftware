#include "PdbMapIntVecFloat.hh"
#include <phool/phool.h>

#include <iostream>

using namespace std;

int verbosity = 0;

PdbMapIntVecFloat::~PdbMapIntVecFloat()
{
  TheMap.clear();
  return;
}

void 
PdbMapIntVecFloat::get_vector(int key, std::vector<float> &fvec) const
{
  map<int, vector<float> >::const_iterator miter = TheMap.find(key);
  fvec.clear();
  if (miter != TheMap.end())
    {
      fvec = miter->second;
    }
  else
    {
      cout << PHWHERE << "No entry for key value " << key << endl;
    }
  return;
}

int
PdbMapIntVecFloat::add_vector(int key, std::vector<float> &fvec)
{
  int iret = 0;
  map<int, vector<float> >::iterator miter = TheMap.find(key);
  if (miter != TheMap.end())
    {
      cout << PHWHERE << "Overwriting entry for key value " << key << endl;
      TheMap.erase(miter);
      iret = 1;
    }
  TheMap[key] = fvec;
  return iret;
}

void
PdbMapIntVecFloat::Verbosity(const int i)
{
  verbosity = i;
  return;
}

void PdbMapIntVecFloat::Clear(Option_t *) { TheMap.clear(); }

void 
PdbMapIntVecFloat::print() const
{
  if (TheMap.empty())
    {
      cout << "No Entries in Map" << endl;
      return;
    }
  map<int, vector<float> >::const_iterator miter;
  vector<float>::const_iterator viter;
  for (miter = TheMap.begin(); miter != TheMap.end(); ++miter)
    {
      cout << "key: " << miter->first << endl;
      int i=0;
      for (viter = (miter->second).begin(); viter != (miter->second).end(); ++viter)
	{
	  cout << "entry " << i << " val: " << *viter << endl;
	  i++;
	}
      cout << "----------" << endl << endl;
    }
  return;
}
