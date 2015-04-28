#include <FROG.h>
#include <sys/stat.h>
#include <string>
#include <pgsearch.h>
#include <dCachesearch.h>

#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

extern char **environ;
static const char sep = ':';
 

const char *
FROG::location(const char * logical_name)
{
  const char * notfound = "";
  int n = 0, i = 0, envset = 0;
  string en, en1, en2;

  if (strcmp(logical_name,"") == 0)
    {
      return notfound;
    }
  
  if(strncmp(logical_name,"/",1) == 0) 	 
    { 	 
      return logical_name; 	 
    }
  
  while (envset == 0 && environ[i])
    {
      en = environ[i];
      if (en.substr(0, 12) == "GSEARCHPATH=")
	{
	  envset = 1;
	  en1 = en.substr(12,en.length());
	  if(en1.empty())
	    {
	      cout << "GSEARCHPATH is an empty string" << endl;
	      exit(1);
	    }
	  while (!en1.empty())
	    {
	      n = en1.find_first_of(sep);
	      if (n != -1)
		{
		  en2 = en1.substr(0,n);
		  en1 = en1.substr(n+1,en1.length());
		}
	      else
		{
		  en2 = en1;
		  en1 = "";
		}
	      if(en2.substr(0,4) == "OBJY")
		{
		  cout << "Objy search is deprecated, please remove OBJY from your GSEARCHPATH env" << endl;
		}
	      else if(en2.substr(0,2) == "PG")
		{
                  pfn.erase(); // erase string from previous calls
		  searchPG(logical_name, pfn);
		  if (!pfn.empty())
		    {
		      return pfn.c_str();
		    }
		}
	      else if(en2.substr(0,6) == "DCACHE")
		{
                  pfn.erase(); // erase string from previous calls
		  searchDC(logical_name, pfn);
		  if (!pfn.empty())
		    {
		      return pfn.c_str();
		    }
		}
	      else 
		{
		  string tem = en2+"/"+logical_name;
		  pfn = localSearch(tem);
		  if (!pfn.empty())
		    {
		      return pfn.c_str();
		    }
		}
	    }
	}
      i++;
    }

  return logical_name; 

}

const char *
FROG::localSearch(const string &logical_name)
{
  struct stat64 stbuf;
  const char * found = "";
  if (stat64(logical_name.c_str(), &stbuf) != -1)
    {
      if ((stbuf.st_mode & S_IFMT) == S_IFREG)
        {
          return logical_name.c_str();
        }
    }
  return found;
}


void
FROG::searchPG(const char * logical_name, string &cp)
{ 
  // Now pgsearches are already known object types...
  pgsearch PG;
  string lname = logical_name;
  //  string cpn = cp;
  PG.search(lname, cp);
  //   cout << "strlen cp: " << strlen(cp) 
  //        << ", sizeof(cp) " << sizeof(cp) << endl;
  //   sprintf(cp,"%s",cpn.c_str());
}

void
FROG::searchDC(const char * logical_name, string &cp)
{ 
  dCachesearch dC;
  string lname = logical_name;
  dC.search(lname, cp);
}
