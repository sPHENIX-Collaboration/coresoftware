///////////////////////////////////////////////////////////////
//
// SpinDBOutput class
// Author      : D. Loomis, D. Neff (from Y. Fukao PHENIX class)
// Description : Utility to read data from spin database
// Created     : 2024-05-12
//
//
///////////////////////////////////////////////////////////////

#ifndef USPIN_SPINDBOUTPUT_H
#define USPIN_SPINDBOUTPUT_H

#include <map>
#include <memory>  // for unique_ptr
#include <string>
#include <vector>

namespace odbc
{
  class Connection;
  class Statement;
  class ResultSet;
};  // namespace odbc

class SpinDBContent;

class SpinDBOutput
{
 public:
  SpinDBOutput() { Initialize(); }
  SpinDBOutput(const char *user)
  {
    Initialize();
    SetUserName(user);
  }
  virtual ~SpinDBOutput() = default;
  void Initialize();
  void SetUserName(const std::string &user)
  {
    user_name = user;
    return;
  }
  void SetDBName(const std::string &dbname);
  void SetTableName(const std::string &tname);
  int PrintDBColumn();
  int PrintDBRawContent(int runnum);
  int PrintDBRawContent(int runnum, int qa_level);
  int CheckRunRow(int runnum);
  int CheckRunRow(int runnum, int qa_level);
  int CheckRunRowStore(int runnum);
  int StoreDBContent(int run1, int run2);
  int StoreDBContent(int run1, int run2, int qa_level);
  void ClearDBContent();
  int GetDBContent(SpinDBContent *&spin_cont, int runnum);
  int GetDBContent(SpinDBContent *&spin_cont, int runnum, int qa_level);
  int GetDBContentStore(SpinDBContent *&spin_cont, int runnum);
  static int CopyDBContent(SpinDBContent &spin_cont, SpinDBContent &spin_cont_copy);
  int GetDefaultQA(int runnum);
  void Verbosity(int verbose = 0) { verbosity = verbose; }

 private:
  static constexpr int ERROR_VALUE{-999};

  int verbosity{0};

  std::string db_name;
  std::string user_name;
  std::string table_name;

  std::map<int, std::unique_ptr<SpinDBContent>> spin_cont_store;

  odbc::Connection *ConnectDB(void);
  int GetDBContent(SpinDBContent &spin_cont, odbc::ResultSet *rs);
  int GetArray(odbc::ResultSet *rs, const std::string &name, std::vector<std::string> &value) const;
  int GetArray(odbc::ResultSet *rs, const std::string &name, float *value, int nvalue);
  int GetArray(odbc::ResultSet *rs, const std::string &name, unsigned int *value, int nvalue);
  int GetArray(odbc::ResultSet *rs, const std::string &name, int *value, int nvalue);
  int GetArray(odbc::ResultSet *rs, const std::string &name, long long *value, int nvalue);
};

#endif /* USPIN_SPINDBOUTPUT_H */
