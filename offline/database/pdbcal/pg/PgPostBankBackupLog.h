// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: PgPostBankBackupLog.hh,v 1.2 2014/05/19 17:06:23 jinhuang Exp $

/*!
 * \file PgPostBankBackupLog.hh
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2014/05/19 17:06:23 $
 */

#ifndef PDBCALPG_PGPOSTBANKBACKUPLOG_H
#define PDBCALPG_PGPOSTBANKBACKUPLOG_H

#include <string>

class TSQLConnection;
class TSQLPreparedStatement;

/*!
 * \brief PgPostBankBackupLog
 */
class PgPostBankBackupLog
{
 public:
  PgPostBankBackupLog(const std::string& TableName, const std::string& Tag);
  virtual ~PgPostBankBackupLog();

 public:
  //! Sets the verbosity of this module (0 by default=quiet).
  virtual void
  Verbosity(const int ival)
  {
    verbosity = ival;
  }

  //! Gets the verbosity of this module.
  virtual int
  Verbosity() const
  {
    return verbosity;
  }

  enum enu_ops
  {
    //! kOptFailed can be multiplied to other operations to flag a failure
    kOptFailed = -1,
    //! positive value means fine
    kOptSuccess = 1,

    kOptNull = 0,

    //! use in PgPostBankBackupManager::fetchAllBank2TFile
    kOptBackup2File = 10,

    //! use in PgPostBankBackupManager::commitAllBankfromTFile
    kOptFile2Db = 20,
    kOptFile2Db_Skip = 21,

    //! use in PgPostBankBackupManager::dumpTable
    kOptDump = 30,

    //! use in PgPostBankBackupManager::CleanTable
    kOptDelete = 100
  };

  //! initialization. Always called to check connections in Log()
  void
  Init();

  //! enter one log entry
  void
  Log(const int rid, const enu_ops ops);

 protected:
  //! The verbosity level. 0 means not verbose at all.
  int verbosity;

  typedef TSQLConnection* TSQLConnection_PTR;
  static TSQLConnection_PTR con;

  std::string tablename;
  std::string tag;

  TSQLPreparedStatement* pstmt;
};

#endif /* PDBCAL_PG_PGPOSTBANKBACKUPLOG_H */
