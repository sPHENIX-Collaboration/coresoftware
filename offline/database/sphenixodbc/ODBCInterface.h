#ifndef ODBC_ODBCINTERFACE_H
#define ODBC_ODBCINTERFACE_H

#include <cstdint>
#include <map>
#include <string>

namespace odbc
{
  class Connection;
  class Statement;
}  // namespace odbc

class ODBCInterface
{
 public:
  static ODBCInterface *instance();

  virtual ~ODBCInterface();

  void Verbosity(uint64_t i) {m_Verbosity = i;}
  uint64_t Verbosity() const {return  m_Verbosity;}
  
  void Print(const std::string & /*what*/ = "ALL") const;

  odbc::Connection *getDBConnection(const std::string &dbname);
  odbc::Statement *getStatement(const std::string &dbname);

  int ConnectionTries() const {return m_ConnectionTries;}
  int SleepMS() const {return m_SleepMS;}

  void Disconnect();
  const std::map<std::string, int>& getNumConnection() const {return m_NumConnection;}
  const std::map<std::string, int>& getNumStatementUse() const {return m_NumStatementUse;}
  
private:

  ODBCInterface() = default;
  static ODBCInterface *__instance;
  int m_ConnectionTries {0};
  int m_SleepMS {0};
  uint64_t m_Verbosity {0};
  static constexpr int m_MAX_NUM_RETRIES = 3000;
  static constexpr int m_MIN_SLEEP_DUR = 200;   // milliseconds
  static constexpr int m_MAX_SLEEP_DUR = 3000;  // milliseconds

  std::map<std::string, odbc::Connection *> m_OdbcConnectionMap;
  std::map<std::string, odbc::Statement *> m_OdbcStatementMap;
  std::map<std::string, int> m_NumConnection;
  std::map<std::string, int> m_NumStatementUse;
};

#endif
