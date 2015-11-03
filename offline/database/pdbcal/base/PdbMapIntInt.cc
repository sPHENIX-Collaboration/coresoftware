#include "PdbMapIntInt.hh"
#include <phool/phool.h>

#include <iostream>

using namespace std;


PdbMapIntInt::PdbMapIntInt():
  verbosity(0)
 {}

PdbMapIntInt::~PdbMapIntInt() { TheMap.clear(); }

void PdbMapIntInt::set_map(std::map<int,int> *inmap) {

  if(verbosity) std::cout << "PdbMapIntInt::set_map() input size = " << inmap->size() << std::endl;
  int count=0;
  map<int,int>::const_iterator miter;
  for (miter = inmap->begin(); miter != inmap->end(); ++miter) {
    //std::cout << count << " " << miter->first << " " << miter->second << std::endl;
    int tmp1 = miter->first;
    int tmp2 = miter->second;
    //std::cout << "PdbMapIntInt::set_map() filling " << tmp1 << " " << tmp2 << std::endl;
    //std::cout << "PdbMapIntInt::set_map() current map size = " << TheMap.size() << std::endl;
    TheMap[tmp1] = tmp2;
    //std::cout << "PdbMapIntInt::set_map() current map size(2) = " << TheMap.size() << std::endl;
    count++;
  }
  if(verbosity) std::cout << "PdbMapIntInt::set_map() output size = " << TheMap.size() << std::endl;

  return;
}

void PdbMapIntInt::Clear(Option_t *) { TheMap.clear(); }

void PdbMapIntInt::print() const {
  if (TheMap.empty())
    {
      cout << "No Entries in Map" << endl;
      return;
    }
  else 
    {
      cout << "Map Size = " << TheMap.size() << endl;
    }
  return;
}


