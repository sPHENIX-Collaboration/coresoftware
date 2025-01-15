// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FROG_FROG_H
#define FROG_FROG_H

#include <map>
#include <string>

namespace odbc
{
  class Connection;
}

class FROG
{
 public:
  FROG() = default;
  virtual ~FROG() = default;

  const char *location(const std::string &logical_name);
  bool localSearch(const std::string &lname);
  bool dCacheSearch(const std::string &lname);
  bool XRootDSearch(const std::string &lname);
  bool LustreSearch(const std::string &lname);
  bool MinIOSearch(const std::string &lname);
  bool RawDataSearch(const std::string &lname);
  bool HpssRawDataSearch(const std::string &lname);
  bool PGSearch(const std::string &lname);
  void Verbosity(const int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  odbc::Connection *GetConnection(const std::string &database);
  void Disconnect();
  static const int m_MAX_NUM_RETRIES{3000};
  static const int m_MIN_SLEEP_DUR{5000};   // milliseconds
  static const int m_MAX_SLEEP_DUR{30000};  // milliseconds

  std::map<std::string, odbc::Connection *> m_OdbcConnectionMap;
  int m_Verbosity{0};
  std::string pfn;
};

#endif
