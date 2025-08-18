//////////////////////////////////////////////////////////////
//
// SpinDBInput class
// Author      : D. Loomis, D. Neff (from Y. Fukao PHENIX class)
// Description : Utility to write data into the spin database
// Created     : 2024-05-12
//
//////////////////////////////////////////////////////////////

#ifndef USPIN_SPINDBINPUT_H
#define USPIN_SPINDBINPUT_H

#include <stdio.h>
#include <sstream>

namespace odbc
{
  class Connection;
  class Statement;
  class ResultSet;
};  // namespace odbc

class SpinDBContent;

class SpinDBInput
{
 public:
  SpinDBInput() { Initialize(); }
  virtual ~SpinDBInput();
  void Initialize();
  int IsConnected();
  int CheckRunRow(int runnum, int qa_level, const std::string &opt = "");
  int CreateRunRow(int runnum, int qa_level);
  int DeleteRunRow(int runnum, int qa_level);
  int InitializeRunRow(SpinDBContent &spin_cont);
  int UpdateDBContent(SpinDBContent &spin_cont);
  int SetDefaultQA(SpinDBContent &spin_cont);
  int UpdateValue(int runnum, int qa_level, const std::string &name, int value);
  int UpdateValue(int runnum, int qa_level, const std::string &name, float value);
  int UpdateValue(int runnum, int qa_level, const std::string &name, double value);
  int UpdateArray(int runnum, int qa_level, const std::string &name, int *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const std::string &name, float *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const std::string &name, double *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const std::string &name, unsigned int *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const std::string &name, long long *value, int nvalue);
  int InitializeValue(int runnum, int qa_level, const std::string &name);
  int InitializeArray(int runnum, int qa_level, const std::string &name, int nvalue);

 private:
  static constexpr std::string DB_NAME{"spinDB_write"};
  static constexpr std::string TABLE_NAME{"spin"};
  static constexpr int ERROR_VALUE{-999};

  odbc::Connection *con{nullptr};
  int run_check;
  int qa_check;

  int UpdateValue(int runnum, int qa_level, const std::string &cmd);
  template <class T>
  int UpdateValueTemp(int runnum, int qa_level, const std::string &name, T value);
  template <class T>
  int UpdateArrayTemp(int runnum, int qa_level, const std::string &name, T *value, int nvalue);
};

#endif /* USPIN_SPINDBINPUT_H */
