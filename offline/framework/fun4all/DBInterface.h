#ifndef FUN4ALL_DBINTERFACE_H
#define FUN4ALL_DBINTERFACE_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class ODBCInterface;
namespace odbc
{
  class Connection;
  class Statement;
}  // namespace odbc

class DBInterface : public SubsysReco
{
 public:
  static DBInterface *instance();

  ~DBInterface() override;

  int process_event(PHCompositeNode *) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *) override;

  void Print(const std::string & /*what*/ = "ALL") const override;

  odbc::Connection *getDBConnection(const std::string &dbname);
  odbc::Statement *getStatement(const std::string &dbname);
  
 private:

  DBInterface(const std::string &name = "DBInterface");
  static DBInterface *__instance;
  ODBCInterface *m_ODBC{nullptr};
};

#endif
