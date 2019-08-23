// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: PgPostBankBackupManager.hh,v 1.8 2014/05/19 17:06:23 jinhuang Exp $

/*!
 * \file PgPostBankBackupManager.hh
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.8 $
 * \date $Date: 2014/05/19 17:06:23 $
 */

#ifndef PDBCALPG_PGPOSTBANKBACKUPMANAGER_H
#define PDBCALPG_PGPOSTBANKBACKUPMANAGER_H

#include <phool/PHTimeStamp.h>

#include <iostream>
#include <string>
#include <vector>

class PdbCalChan;
class PgPostBankBackupStorage;
class TSQLStatement;
class TSQLPreparedStatement;
class TSQLResultSet;

/*!
 * \brief PgPostBankBackupManager
 */
class PgPostBankBackupManager
{
 public:
  PgPostBankBackupManager(const std::string &Tag = "");
  virtual ~PgPostBankBackupManager();

  //! list of int
  typedef std::vector<int> rid_list_t;

 public:
  //! check whether an record with certain rid already there
  bool
  isRIdExist(const std::string &bankName, int rid);

  //! total row counts
  int getTotalRowCount(const std::string &bankName);

  //! list rid in a bank
  rid_list_t
  getListOfRId(const std::string &bankName,
               const std::string &condition = "");

  //! Build PgPostBankBackupStorage from one database entry
  //! using code should delete the returned object
  //! \param[in] bankName can be bank name (e.g. calib.fvtx.alignment) or pSQL table name (e.g. calibfvtxalignment)
  PgPostBankBackupStorage *
  fetchBank(const std::string &bankName, int rid);

  //! Use PgPostBankBackupStorage to restore the original entry in database
  //! \param[in] bs PgPostBankBackupStorage object. This function will not own this object
  //! \return true if success
  bool
  commit(PgPostBankBackupStorage *bs);

  //! TFile -> Database for official operation use
  int commitAllBankfromTFile(const std::string &in_file);

  //! Database -> TFile for official operation use
  int fetchBank2TFile(const std::string &bankName, const rid_list_t &rid_list,
                      const std::string &out_file);

  //! Database -> TFile for all all records satisifying  record_selection
  int fetchAllBank2TFile(const std::string &bankName,
                         const std::string &record_selection, const std::string &out_file_base,
                         int splitting_limit = 100000);

  //! dump databse to text
  void
  dumpTable(const std::string &bankName, std::ostream &out);

  //! remove records which is NOT used for reconstruction.
  //! Records with insert timeafter min_save_time will be preserved
  int CleanTable(const std::string &bankName,
                 const PHTimeStamp &min_preserve_time = PHTimeStamp::PHFarFuture,
                 const bool do_delete = false, const std::string &log_file_name = "");

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

  static std::string
  getBankBaseName(const std::string &bank_classname);

  //! clear list of TSQLStatement
  static void
  deleteSQLStatement(TSQLStatement *stmt);

  //! clear list of TSQLStatement
  static void
  deleteODBCPreparedStatement(TSQLPreparedStatement *stmt);

  //! Hash the object
  static ULong_t
  HashPdbCalChan(const PdbCalChan &c);

 protected:
  PgPostBankBackupStorage *
  SQLResultSet2BackupStorage(TSQLResultSet *rs, const std::string &table_name);

  //! The verbosity level. 0 means not verbose at all.
  int verbosity;

  std::string tag;
};

#endif /* PDBCAL_PG_PGPOSTBANKBACKUPMANAGER_H */
