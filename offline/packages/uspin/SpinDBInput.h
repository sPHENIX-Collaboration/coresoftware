//////////////////////////////////////////////////////////////
//
// SpinDBInput class
// Author      : D. Loomis (from Y. Fukao PHENIX class)
// Description : Utility to write data into the spin database
// Created     : 2024-05-12
//
//////////////////////////////////////////////////////////////

#ifndef USPIN_SPINDBINPUT_H
#define USPIN_SPINDBINPUT_H

#include "SpinDBContent.h"

#include <stdio.h>
#include <sstream>

namespace odbc
{
  class Connection;
  class Statement;
  class ResultSet;
};  // namespace odbc

class SpinDBInput
{
 public:
  SpinDBInput() { Initialize(); }
  virtual ~SpinDBInput();
  void Initialize();
  int IsConnected();
  int CheckRunRow(int runnum, int qa_level, const char *opt = "");
  int CreateRunRow(int runnum, int qa_level);
  int DeleteRunRow(int runnum, int qa_level);
  // int CheckQARunRow(int runnum);
  // int CreateQARunRow(int runnum);
  // int DeleteQARunRow(int runnum);
  // int SetQADefault(int runnum,int qa_level);
  int InitializeRunRow(SpinDBContent spin_cont);
  int UpdateDBContent(SpinDBContent spin_cont);
  int UpdateValue(int runnum, int qa_level, const char *name, int value);
  int UpdateValue(int runnum, int qa_level, const char *name, float value);
  int UpdateValue(int runnum, int qa_level, const char *name, double value);
  int UpdateArray(int runnum, int qa_level, const char *name, int *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const char *name, float *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const char *name, double *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const char *name, unsigned int *value, int nvalue);
  int UpdateArray(int runnum, int qa_level, const char *name, long long *value, int nvalue);
  int InitializeValue(int runnum, int qa_level, const char *name);
  int InitializeArray(int runnum, int qa_level, const char *name, int nvalue);

 private:
  static const std::string DB_NAME;
  static const std::string TABLE_NAME;
  static const int ERROR_VALUE;

  odbc::Connection *con{nullptr};
  int run_check;
  int qa_check;

  int UpdateValue(int runnum, int qa_level, const char *cmd);
  template <class T>
  int UpdateValueTemp(int runnum, int qa_level, const char *name, T value);
  template <class T>
  int UpdateArrayTemp(int runnum, int qa_level, const char *name, T *value, int nvalue);
};

#endif /* USPIN_SPINDBINPUT_H */
