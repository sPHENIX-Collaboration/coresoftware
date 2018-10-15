// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FROG_FROG_H
#define FROG_FROG_H

#include <string>
namespace odbc
{
 class Connection;
}

class FROG
{
 public:
  FROG();
  virtual ~FROG() {}

  const char *location(const std::string &logical_name);
  bool localSearch(const std::string &lname);
  bool dCacheSearch(const std::string &lname);
  bool PGSearch(const std::string &lname);
  void Verbosity(const int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  bool GetConnection();
  void Disconnect();
  std::string pfn;
  int m_Verbosity;
  odbc::Connection *m_OdbcConnection;
};

#endif
