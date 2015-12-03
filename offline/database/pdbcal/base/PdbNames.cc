// This class parses the BankName and splits it up into pieces that
// are used to determine which directory the object should be placed,
// the name of the database, and the name of the container inside the
// database.  It should evolve to support A naming scheme where the
// database is given an additional four digit number "0000" that will
// allow for the database to be broken up if the size exceeds a
// predetermined threshold.
							       
// Author: Mark Pollack                                        

#include <cstdio>
#include <string>
#include <sstream>
#include "PdbNames.h"
#include "phool/PHPointerList.h"
#include "PHString.h"

using namespace std;

PdbNames::PdbNames(const PHString& newBankName, 
		   const PHString& bootFile, 
		   int dbFileNumber)
{
  bootFileFullPath = bootFile;
  setBankName(newBankName, dbFileNumber);
}

PHString 
PdbNames::getDbFileName(const PHString& sysName) const
{
  PHPointerList<PHString> directories;
  bootFileFullPath.split(directories, "/");
  PHString dbFileFullPath;
  for (size_t i = 0; i < directories.length() - 1; i++)
    {
      if (directories[i]->length())
	{
	  if ( ! directories[i]->find("::") )
	    {
	      dbFileFullPath = dbFileFullPath + "/" + *(directories[i]);
	    }
	}
    }
  dbFileFullPath = dbFileFullPath + "/" + sysName + ".pdb";
  
  directories.clearAndDestroy();
  return dbFileFullPath;
}

void 
PdbNames::setBankName(const PHString & newBankName, 
		      int dbFileNumber)
{
  bankName = newBankName;
  PHPointerList<PHString> names;
  bankName.split(names, ".");
  if (names.length() >= 3) 
    {
      calType       = *(names[0]);
      detType       = *(names[1]);
      char filebuf[5];
      sprintf(filebuf, "%04d", dbFileNumber);
      tagDbSysName  = bankName + "Tag" + filebuf;
      calDbSysName  = bankName + "Cal" + filebuf;
      tagDbFileName = getDbFileName(tagDbSysName);
      calDbFileName = getDbFileName(calDbSysName);
      tagContName   = detType + "TagContainer";
      calContName   = detType + "CalContainer";
    }
  else 
    {
      cout << "PdbNames::setBankName(...)" << endl
	   << "\tError" << endl
	   << "\tbank name '" << newBankName << "' could not be split" << endl;
    }
  names.clearAndDestroy();
}

void
PdbNames::incrementCalContainerName()
{
   std::string s(calContName.getString());
   std::string c("CalContainer");
   std::stringstream istr; 
   std::string istring;
   int i = 1;

   std::string::size_type n;
   if ((n = s.rfind(c)) == std::string::npos)
     {
       return;
     }
   if (s.length() > n + c.length())
     {
       std::stringstream tmp(s.substr(n + c.length(), s.length()));
       tmp >> i;
       s.erase(n + c.length(), s.length());
       i++;
     }

   istr << i;
   istring = istr.str();
   s.insert((n + c.length()), istring);
   calContName = s.c_str();

}
