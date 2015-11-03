//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 2004
//
//  Implementation of class PdbParameter
//
//  Author: Cesar

#include <PdbParameterError.hh>

#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

PdbParameterError::PdbParameterError()
{
  thePar = 0;
  theParError = 0;
  setName("UNKNOWN");
}

PdbParameterError::PdbParameterError(const float value)
{
  thePar = value;
  theParError = 0;
  setName("UNKNOWN");
}

PdbParameterError::PdbParameterError(const float value, const float error)
{
  thePar = value;
  theParError = error;
  setName("UNKNOWN");
}

PdbParameterError::PdbParameterError(const float value, const float error, const char *name)
{
  thePar = value;
  theParError = error;
  setName(name);
}

void
PdbParameterError::print() const
{
  cout << theName << ": " << thePar << " +/- " << theParError << endl;
}

void 
PdbParameterError::setName(const char *name)
{
  // this construct protects against not zero terminated strings
  // it copies only a fixed number of chars into the theName array
  // (sizeof(theName)-1 leaves space for a zero terminated char
  // at the end. strlen has problems when used with non zero
  // terminated strings and might crash
  // this code might still crash in the cout of name, but that is 
  // already in the exit() part
  // Flawfinder: ignore signals flawfinder to ignore this line, strncpy does not zero terminate the string
  // if the length is exceeded, the termination is done in the line afterwards where the last
  // character is set to \0
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
