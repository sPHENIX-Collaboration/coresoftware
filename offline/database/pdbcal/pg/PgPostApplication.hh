#ifndef PGPOSTAPPLICATION_H
#define PGPOSTAPPLICATION_H

#include <pdbcalbase/PdbApplication.hh>
#include <phool/phool.h>
#include <pdbcalbase/Pdb.hh>

#include <string>

class TSQLConnection;

class PgPostApplication : public PdbApplication{ 

 protected: 
  PgPostApplication(const std::string &dbname); 

 public: 
  virtual ~PgPostApplication(); 
  PdbStatus startUpdate();
  PdbStatus startRead();
  PdbStatus commit();
  PdbStatus commit(PdbCalBank *);
  PdbStatus commit(PdbCalBank *, int rid, long,long,long);

  PdbStatus abort();
  PdbStatus isActive() { return 0; }
  TSQLConnection * getConnection();
  size_t getTagFileSize(const char *)  { return 0; }
  size_t getCalFileSize(const char * ) { return 0; }

  static int Register(const std::string &dbname = "calibrations");
  static int releaseConnection( );
  static PgPostApplication *instance();
  int setDBName(const char *name);
  int DisconnectDB();

 protected:
  static TSQLConnection *        con;
  static PgPostApplication *mySpecificCopy;
  bool readOnly;
  std::string dsn;
}; 



#endif /* PGPOSTAPPLICATION_H */ 


