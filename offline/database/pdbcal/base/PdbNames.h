#ifndef __PDBNAMES_HH__
#define __PDBNAMES_HH__

#include "PHString.h"

class PdbNames {
public:
  PdbNames() {}
   PdbNames(const PHString& bankName, 
	    const PHString& bootFile, 
	    int dbFileNumber = 0);
 
  virtual ~PdbNames() {}

   const char* getCalibrationType()  const { return calType.getString(); }
   const char* getDetectorType()     const { return detType.getString(); }
   const char* getTagDbSysName()     const { return tagDbSysName.getString(); }
   const char* getCalDbSysName()     const { return calDbSysName.getString(); }
   const char* getTagContainerName() const { return tagContName.getString(); }
   const char* getCalContainerName() const { return calContName.getString(); }
   const char* getTagDbFileName()    const { return tagDbFileName.getString(); }
   const char* getCalDbFileName()    const { return calDbFileName.getString(); }
   
  PHString getDbFileName(const PHString &) const;
  void incrementCalContainerName();

private:
  void setBankName(const PHString& bankName,int dbFileNumber);
   
private:
   PHString bankName;   
   PHString calType; 
   PHString detType; 
   PHString tagDbSysName;  
   PHString calDbSysName;  
   PHString tagDbFileName;
   PHString calDbFileName;
   PHString tagContName;
   PHString calContName;
   PHString bootFileFullPath;
};

#endif /* __PDBNAMES_HH__ */
