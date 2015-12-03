//  Implementation of class PdbIndex
//  Author: Federica Messer

#include <PdbIndex.h>
#include <PHString.h>

#include <cstdlib>
#include <cstring>
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

void 
PdbIndex::setName(const char name[])
{
  // this construct protects against not zero terminated strings
  // it copies only a fixed number of chars into the theName array
  // (sizeof(theName)-1 leaves space for a zero terminated char
  // at the end. strlen has problems when used with non zero
  // terminated strings and might crash
  // this code might still crash in the cout of name, but that is 
  // already in the exit() part
  strncpy(theName,name,(sizeof(theName)-1)); /* Flawfinder: ignore */
  theName[sizeof(theName)-1] = '\0';
  if (strncmp(theName,name,sizeof(theName)))
    {
      cout << "Name exceeds maximum length of "
	   << (sizeof(theName)-1) 
           << " characters or is not zero terminated" << endl;
      cout << "Max length name: " << theName << endl;
      cout << "There is no point in continuing, fix your code and try again" << endl;
      cout << "Name used (code might crash now when printing out not zero terminated string): " << name << endl;
	exit(1);
    }
}

PHBoolean 
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
