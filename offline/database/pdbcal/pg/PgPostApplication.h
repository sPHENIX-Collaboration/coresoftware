// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PDBCALPG_PGPOSTAPPLICATION_H
#define PDBCALPG_PGPOSTAPPLICATION_H

#include <pdbcalbase/Pdb.h>
#include <pdbcalbase/PdbApplication.h>

#include <string>

class PdbCalBank;
class TSQLConnection;

class PgPostApplication : public PdbApplication
{
 protected:
  PgPostApplication(const std::string &dbname);

 public:
  virtual ~PgPostApplication();
  PdbStatus startUpdate();
  PdbStatus startRead();
  PdbStatus commit();
  PdbStatus commit(PdbCalBank *);
  PdbStatus commit(PdbCalBank *, int rid, long, long, long);

  PdbStatus abort();
  PdbStatus isActive() { return 0; }
  TSQLConnection *getConnection();

  static int Register(const std::string &dbname = "calibrations");
  static int releaseConnection();
  static PgPostApplication *instance();
  int setDBName(const std::string &name);
  int DisconnectDB();

 protected:
  static TSQLConnection *con;
  static PgPostApplication *mySpecificCopy;
  bool readOnly;
  std::string dsn;
};

#endif /* PDBCAL_PG_PGPOSTAPPLICATION_H */
