#include "PgPostApplication.h"
#include "PgPostBankWrapper.h"
#include "PgPostBankWrapperManager.h"

#include <pdbcalbase/PdbApplicationFactory.h>
#include <pdbcalbase/PHGenericFactoryT.h>
#include <pdbcalbase/Pdb.h>
#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbCalBank.h>

#include <phool/phool.h>

#include <RDBC/TSQLDriverManager.h>
#include <RDBC/TSQLConnection.h>

#include <iostream>
#include <memory>
#include <sstream>

using namespace std;

PgPostApplication *PgPostApplication::mySpecificCopy = nullptr;
TSQLConnection *PgPostApplication::con = nullptr;

namespace
{
PdbApplication *singletonCreator()
{
  // rememeber that this will not neccessarily return a
  // pointer to the singleton PgPostApplication. If
  // an Objy application is instantiated, it will return 0.

  cout << "*** deprecated interface -- use PdbApplication::instance()" << endl;
  return PgPostApplication::instance();
}

const std::string name = "Pg";
const bool registered =
    PdbApplicationFactory::instance().registerCreator(name, singletonCreator, "PdbApplication");
}  // namespace

PgPostApplication *PgPostApplication::instance()
{
  return mySpecificCopy;
}

int PgPostApplication::Register(const string &dbname)
{
  if (__instance.get()) return -1;
  mySpecificCopy = new PgPostApplication(dbname);
  __instance = std::unique_ptr<PdbApplication>(mySpecificCopy);
  return 0;
}

int PgPostApplication::releaseConnection()
{
  if (!__instance.get()) return -1;
  if (con)
  {
    con->Close();
    con = nullptr;
  }
  return 0;
}

PgPostApplication::PgPostApplication(const string &dbname)
  : readOnly(1)
  , dsn(dbname)
{
}

PgPostApplication::~PgPostApplication()
{
  mySpecificCopy = 0;
  if (con)
  {
    con->Close();
    con = 0;
  }
}

TSQLConnection *PgPostApplication::getConnection()
{
  if (!__instance.get())
  {
    PgPostApplication::Register();
  }
  if (!con)
  {
    ostringstream constr;
    constr << "dsn=" << dsn << "; uid=phnxrc; pwd= ";
    con = gSQLDriverManager->GetConnection(constr.str().c_str());
  }
  return con;
}

int PgPostApplication::setDBName(const string &name)
{
  if (dsn != name)
  {
    dsn = name;
    if (con)
    {
      con->Close();
    }
    ostringstream constr;
    constr << "dsn=" << dsn << "; uid=phnxrc; pwd= ";
    con = gSQLDriverManager->GetConnection(constr.str().c_str());
    if (!con)
    {
      cout << PHWHERE << " Could not open data base " << dsn << endl;
    }
  }
  return 0;
}

PdbStatus
PgPostApplication::startUpdate()
{
  readOnly = 0;
  return 1;
}

PdbStatus
PgPostApplication::startRead()
{
  readOnly = 1;
  return 1;
}

PdbStatus
PgPostApplication::abort()
{
  PgPostBankWrapperManager::instance().clear();
  // will roll back a Xact in progres
  if (con) con->Rollback();
  return 1;
}

PdbStatus
PgPostApplication::commit()
{
  if (readOnly)
  {
    // Forget about the registered wrappers.
    PgPostBankWrapperManager::instance().clear();
    return 1;
  }
  else
  {
    // Commit the registered wrappers (this will also forget them).
    bool ok = PgPostBankWrapperManager::instance().commit();
    if (con) con->Commit();
    return (ok == true);
  }
}

PdbStatus
PgPostApplication::commit(PdbCalBank *b)
{
  PgPostBankWrapper *tb = dynamic_cast<PgPostBankWrapper *>(b);
  if (tb)
  {
    return tb->commit();
  }
  else
  {
    cerr << PHWHERE << "Can only commit PgPostBankWrapper2" << endl;
    return 0;
  }
}

PdbStatus
PgPostApplication::commit(PdbCalBank *b, int /*rid*/, long /*it*/, long /*st*/, long /*et*/)
{
  PgPostBankWrapper *tb = dynamic_cast<PgPostBankWrapper *>(b);
  if (tb)
  {
    return tb->commit();
  }
  else
  {
    cerr << PHWHERE << "Can only commit PgPostBankWrapper2" << endl;
    return 0;
  }
}

int PgPostApplication::DisconnectDB()
{
  if (con)
  {
    con->Close();
  }
  con = nullptr;
  return 0;
}
