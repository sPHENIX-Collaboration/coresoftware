//  Implementation of class PdbIndex
//  Author: Federica Messer

#include "PdbIndex.h"

#include <cstdlib>
#include <string>
#include <iostream>

using namespace std;

PdbIndex::PdbIndex()
{
  setMin(0);
  setMax(0);
  setValue(0);
  setName("UNKNOWN");
}

PdbIndex::PdbIndex(const int mini, const int maxi)
{
  setMin(mini);
  setMax(maxi);
  setValue(mini);
  setName("UNKNOWN");
}

PdbIndex::PdbIndex(const int mini, const int maxi, const char name[])
{
  setMin(mini);
  setMax(maxi);
  setValue(mini);
  setName(name);
}

PdbIndex::PdbIndex(const int mini, const int maxi, const int val, const char* name)
{
  setMin(mini);
  setMax(maxi);
  setValue(mini);
  setName(name);
}

void 
PdbIndex::print() const
{
  cout << " PdbIndex Description: \n"
       << "\t Name          :" << theName     << "\n" 
       << "\t Minimus set   :" << theMinimum  << "\n"
       << "\t Maximum set   :" << theMaximum  << "\n"
       << "\t Present Value :" << theValue    << endl;    
}


bool 
PdbIndex::setValue(const int val)
{
  if (val <= theMaximum && val >= theMinimum)
    {
      theValue = val;
      return True;
    }
  else
    {
      theValue = val;
      return False;
    }
}
