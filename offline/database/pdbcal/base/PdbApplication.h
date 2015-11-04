#ifndef __PDBAPPLICATION_H__
#define __PDBAPPLICATION_H__

#include "Pdb.h"

#include <cstddef> // for size_t
#include <memory>  // for auto_ptr

class PdbCalBank;
class PHTimeStamp;

class PdbApplication
{

protected:
  PdbApplication() {}
  virtual ~PdbApplication() {}

public:
  static PdbApplication *instance();


  virtual PdbStatus startUpdate() = 0;
  virtual PdbStatus startRead() = 0;
  virtual PdbStatus commit() = 0;
  virtual PdbStatus abort() = 0;
  virtual PdbStatus isActive() = 0;
  virtual PdbStatus commit(PdbCalBank *) = 0;
  virtual PdbStatus commit(PdbCalBank *,int,long,long,long) = 0;
  //
  // Should return the file (database) size in bytes.
  //
  virtual size_t getTagFileSize(const char *) = 0;
  virtual size_t getCalFileSize(const char *) = 0;

  virtual int setDBName(const char *name) = 0;
  virtual int DisconnectDB() = 0;

protected:
  // Wrap the singleton object in an auto_ptr so it gets cleaned up.
  // Even though auto_ptr is often a horrible choice, here it should be ok
  // since __instance will never be copied. If this wouldn't be visible for
  // CINT a good choice would have been boost::scoped_ptr, or in C++11
  // std::unique_ptr.
  friend class std::auto_ptr<PdbApplication>;
  static std::auto_ptr<PdbApplication> __instance;
};

#endif /* PDBAPPLICATION_H */
