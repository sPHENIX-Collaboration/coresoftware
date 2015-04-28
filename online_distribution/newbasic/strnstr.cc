/* 

This function is modelled after the "strstr" function in
string.h. There is no standard function which locates a substring in a
string which is not 0-terminated. In all other respects the strnstr
function is suuposed to behave as strstr, that is, it returns a
pointer to the first occurence of string s2 in string s1, only that we
can deal here without not 0-terminated strings and rather tell the
length of the strings by the s1len and s2len parameters.  

On some systems the memmem function does this, but it is not available on 
all systems.

*/
#include <iostream>
#include <cstring>

char * strnstr (const char *s1, size_t s1len, const char *s2, size_t s2len)
{

  size_t i;
  char *c;

  /* if s2len is 0, we are done */
  if (s2len <= 0) return (char *) s1;

  /* if s2len > s1len we will not find the substring here. */
  if (s2len > s1len ) return 0;

  char *s2copy = new char[s2len+1];
  c = s2copy;
  for (i=0; i<s2len; i++) *c++ = s2[i];
  *c = 0;

  c = (char *) s1;
  for (i=0; i <= s1len - s2len; i++)
    {
      if (strncmp(c, s2copy, s2len) == 0)
	{
	  delete [] s2copy;
	  return c;
	}
      c++;
    }
  delete [] s2copy;
  return 0;
}
  
