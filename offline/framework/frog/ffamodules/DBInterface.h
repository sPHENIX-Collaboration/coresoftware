#ifndef FFAMODULES_DBINTERFACE_H
#define FFAMODULES_DBINTERFACE_H

#include <map>
#include <string>

namespace odbc
{
  class Connection;
  class Statement;
}  // namespace odbc

class DBInterface
{
 public:
  static DBInterface *instance();

  virtual ~DBInterface();

  odbc::Connection *getDBConnection(const std::string &dbname, int verbosity = 0);
  odbc::Statement *getStatement(const std::string &dbname, int verbosity = 0);
  
 private:

  DBInterface() = default;
  static DBInterface *__instance;
  int m_ConnectionTries {0};
  int m_SleepMS {0};
  static constexpr int m_MAX_NUM_RETRIES = 3000;
  static constexpr int m_MIN_SLEEP_DUR = 200;   // milliseconds
  static constexpr int m_MAX_SLEEP_DUR = 3000;  // milliseconds

  std::map<std::string, odbc::Connection *> m_OdbcConnectionMap;
  std::map<std::string, odbc::Statement *> m_OdbcStatementMap;
  std::map<std::string, int> m_NumConnection;
  std::map<std::string, int> m_NumStatementUse;
};

#endif
