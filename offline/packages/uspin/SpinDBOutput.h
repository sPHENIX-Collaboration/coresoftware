///////////////////////////////////////////////////////////////
//
// SpinDBOutput class
// Author      : D. Loomis (from Y. Fukao PHENIX class)
// Description : Utility to read data from spin database
// Created     : 2024-05-12
//
//
///////////////////////////////////////////////////////////////

#ifndef _SPINDBOUTPUT_
#define _SPINDBOUTPUT_

#include "SpinDBContent.h"

#include <TObject.h>

#include <map>
#include <string>

#define QA_ERROR_VALUE -999

namespace odbc{
  class Connection;
  class Statement;
  class ResultSet;
};

class SpinDBOutput : public TObject{
public:
  SpinDBOutput(void){Initialize();}
  SpinDBOutput(const char *user){Initialize(); SetUserName(user);}
  virtual ~SpinDBOutput(void){Delete();}
  void Initialize(void);
  void SetUserName(const char *user){user_name=user; return;}
  void SetDBName(const char *dbname);
  void SetTableName(const char *tname);
  int PrintDBColumn(void);
  int PrintDBRawContent(int runnum,int qa_level=QA_ERROR_VALUE);
  int CheckRunRow(int runnum,int qa_level=QA_ERROR_VALUE);
  int CheckRunRowStore(int runnum);
  int StoreDBContent(int run1,int run2,int qa_level=QA_ERROR_VALUE);
  void ClearDBContent(void);
  int GetDBContent(SpinDBContent &spin_cont,int runnum,int qa_level=QA_ERROR_VALUE);
  int GetDBContentStore(SpinDBContent &spin_cont,int runnum);
  int GetDefaultQA(int runnum);

private:

  static const int ERROR_VALUE;

  std::string db_name;
  std::string user_name;
  std::string table_name;

  SpinDBContent spin_cont_store1;
  std::map<int,SpinDBContent> spin_cont_store;

  odbc::Connection *ConnectDB(void);
  int GetDBContent(SpinDBContent &spin_cont,odbc::ResultSet *rs);
  int GetArray(odbc::ResultSet *rs,const char *name,std::vector<std::string> &value);
  int GetArray(odbc::ResultSet *rs,const char *name,float *value,int nvalue);
  int GetArray(odbc::ResultSet *rs,const char *name,unsigned int *value,int nvalue);
  int GetArray(odbc::ResultSet *rs,const char *name,int *value,int nvalue);
  int GetArray(odbc::ResultSet *rs,const char *name,long long *value,int nvalue);

  ClassDef(SpinDBOutput,1)
};

#endif /* _SPINDBOUTPUT_ */
