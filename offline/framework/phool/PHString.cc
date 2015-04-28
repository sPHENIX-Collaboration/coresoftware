//  The PHOOL's Software
//  Copyright (C) PHENIX collaboration, 1999
//
//  Implementation of class PHString
//
//  Author: Matthias Messer

#include "PHString.h"
#include "PHPointerList.h"
#include "PHPointerListIterator.h"
#include <cstring>
#include <cstdio>
#include <iostream>

using namespace std;

void PHString::copy_(const char* s)
{
   delete [] string;
   if (!s)
     {
       string = NULL;
     }
   else
     {
       string = new char[strlen(s) + 1];
       strcpy(string, s);
    }
}

PHString::PHString():
  string(NULL)
{
   copy_("");
}

PHString::~PHString()
{
   delete [] string;
}

PHString::PHString(const PHString& phs):
  string(NULL)
{
   copy_(phs.string);
}

PHString::PHString(const char* s):
  string(NULL)
{
   copy_(s);
}

void PHString::setString(char * s)
{
   delete [] string;
   string = s;
}

PHBoolean PHString::find(const PHString& subString) const
{
   if (strstr(string, subString.getString()))
      return True;
   else
      return False;
}

PHBoolean PHString::operator== (const PHString& phs) const
{
   if (strcmp(string, phs.string) == 0)
      return True;
   else
      return False;
}

PHBoolean PHString::operator!= (const PHString& phs) const
{
   if (strcmp(string, phs.string) == 0)
      return False;
   else
      return True;
}

PHString& PHString::operator= (const PHString& phs)
{
   copy_(phs.string);
   return *this;
}

PHString& PHString::operator+= (const PHString& phs)
{
   char *buffer = string;
   string = new char[length()+phs.length()+1];
   if (buffer) {
     strcpy(string, buffer);
   }
   strcat(string, phs.string);
   delete [] buffer;
   return *this;
}

size_t 
PHString::split(PHPointerList<PHString>& list, const char *sep) const
{
  char *buffer = new char[length() + 1];
  char *workStr = buffer;
  char *searchStr = buffer;
  char *ptr = buffer;
  size_t numberOfSubStrings = 0;
  size_t seplen = strlen(sep);
  
  strcpy(workStr, string);
  list.clearAndDestroy();
  while (strlen(searchStr) > 0)
    {
      // Make sure there's enough of the string left to bother searching
      if (strlen(searchStr) <= seplen)
	{
	  break;
	}

      // See if we're looking at an instance of the separator right here.
      if (strncmp(searchStr, sep, seplen) == 0)
	{
	  // Super.  Insert a terminator to mark the end of the
	  // previous bit of string.  This will affect the stuff that
	  // workStr points to too.
	  *searchStr = '\0';
	  list.append(new PHString(workStr));
	  numberOfSubStrings++;
	  workStr = searchStr + seplen;
	  ptr = workStr;
	  searchStr = workStr - 1;
	}
      searchStr++;
    }

  // Some cleaning up.
  if (ptr)
    {
      list.append(new PHString(ptr));
      numberOfSubStrings++;
    }
  
  delete [] buffer;
  return numberOfSubStrings;
}

//
// friend operators and external functions
//
PHString 
operator+ (const PHString& phs1, const PHString& phs2)
{
   PHString newString(phs1);
   newString += phs2;
   return newString;
}

ostream & 
operator << (ostream & s, const PHString & q)
{
  return s << q.getString();
}

PHString 
join(const PHPointerList<PHString>& list, const PHString& seperator)
{
   PHPointerListIterator<PHString> iter(list);
   PHString result = *(iter());
   PHString *element;
   while ((element = iter()))
     {
       result = result + seperator + *element;
     }
   return result;
}
