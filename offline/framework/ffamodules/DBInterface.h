#ifndef FFAMODULES_DBINTERFACE_H
#define FFAMODULES_DBINTERFACE_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

namespace odbc
{
  class  Connection;
}

class DBInterface : public SubsysReco
{
 public:
  static DBInterface *instance();

  ~DBInterface() override = default;

  /// Called at the end of all processing.
  int InitRun(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  double getDVal(const std::string &name);
  static odbc::Connection *getDBConnection(const std::string &dbname, int verbosity = 0);
  int getRunTime(int runnumber);

private:
  union u_value
  {
    uint64_t uidata;
    int64_t idata;
    double ddata;
  u_value(uint64_t in): uidata(in) {}
  u_value(int64_t in): idata(in) {}
  u_value(double in): ddata(in) {}
  };
  
  DBInterface(const std::string &name = "DBInterface");
  static DBInterface *__instance;
  static constexpr int m_MAX_NUM_RETRIES = 3000;
  static constexpr int m_MIN_SLEEP_DUR =  200; // milliseconds
  static constexpr int m_MAX_SLEEP_DUR = 3000; // milliseconds

  std::map<std::string, uint64_t> m_ValueMap;
};

#endif
