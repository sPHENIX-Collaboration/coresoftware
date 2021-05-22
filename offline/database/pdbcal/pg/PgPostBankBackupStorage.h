// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 * \file PgPostBankBackupStorage.hh
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2014/04/14 22:56:48 $
 */

#ifndef PDBCALPG_PGPOSTBANKBACKUPSTORAGE_H
#define PDBCALPG_PGPOSTBANKBACKUPSTORAGE_H

#include <pdbcalbase/PdbCalChan.h>

#include <phool/PHTimeStamp.h>

#include <TNamed.h>
#include <TObject.h>

#include <cstddef>
#include <string>

class PgPostCalBank;
class PdbCalBank;

/*!
 * \brief PgPostBankBackupStorage
 */
class PgPostBankBackupStorage : public TNamed
{
 public:
  //! constructor
  //! \param[in] b pointer to a bank (not wrappers!). PgPostBankBackupStorage will own this pointer
  PgPostBankBackupStorage(PdbCalBank *b);

  PgPostBankBackupStorage();

  virtual ~PgPostBankBackupStorage();

  virtual void
  Print(Option_t *option = "") const;

  bool
  isValid() const;

  void
  format_name_title();

  //! use this storage object to recover the PdbCalBankWrapper/PdbCalBankWrapper2
  //! user code are responsible to delete the returned PdbCalBank
  PgPostCalBank *
  createBank();

  class BankHeader : public TObject
  {
   public:
    enum
    {
      INVALID_BANKID = -999999
    };

    BankHeader()
      : bankID(INVALID_BANKID)
      , description(
            "not assigned")
      , userName("not assigned")
      , tableName(
            "not assigned")
      , rid(INVALID_BANKID)
    {
    }

    ~BankHeader()
    {
    }

    virtual void
    Print(Option_t *option = "") const;

    int
    getBankID() const
    {
      return bankID;
    }

    PHTimeStamp
    getInsertTime() const
    {
      return insertTime;
    }

    PHTimeStamp
    getStartValTime() const
    {
      return startValTime;
    }

    PHTimeStamp
    getEndValTime() const
    {
      return endValTime;
    }

    std::string
    getDescription() const
    {
      return description;
    }

    std::string
    getUserName() const
    {
      return userName;
    }

    std::string
    getTableName() const
    {
      return tableName;
    }

    int
    getRId() const
    {
      return rid;
    }

    void
    setBankID(const int val)
    {
      bankID = val;
    }

    void
    setInsertTime(const PHTimeStamp &val)
    {
      insertTime = val;
    }

    void
    setStartValTime(const PHTimeStamp &val)
    {
      startValTime = val;
    }

    void
    setEndValTime(const PHTimeStamp &val)
    {
      endValTime = val;
    }

    void
    setDescription(const std::string &val)
    {
      description = val;
    }

    void
    setUserName(const std::string &val)
    {
      userName = val;
    }

    void
    setTableName(const std::string &val)
    {
      tableName = val;
    }

    void
    setRId(int r)
    {
      rid = r;
    }

    std::string
    get_id_string() const;

   private:
    int bankID;
    PHTimeStamp insertTime;
    PHTimeStamp startValTime;
    PHTimeStamp endValTime;
    std::string description;
    std::string userName;
    std::string tableName;
    int rid;

    ClassDefOverride(PgPostBankBackupStorage::BankHeader, 1)
  };

  void
  printEntry(size_t s)
  {
    getEntry(s).print();
  }

  size_t
  getLength() const;

  PdbCalChan &
  getEntry(size_t pos);

  const PdbCalBank *
  get_bank() const
  {
    return bank;
  }

  PdbCalBank *
  get_bank()
  {
    return bank;
  }

  const BankHeader &
  get_database_header() const
  {
    return database_header;
  }

  BankHeader &
  get_database_header()
  {
    return database_header;
  }

  void
  set_database_header(const BankHeader &databaseHeader)
  {
    database_header = databaseHeader;
  }

  const std::string &
  get_obj_classname() const
  {
    return obj_classname;
  }

  std::string &
  get_obj_classname()
  {
    return obj_classname;
  }

  void
  set_obj_classname(const std::string &name)
  {
    obj_classname = name;
  }

  const BankHeader &
  get_obj_header() const
  {
    return obj_header;
  }

  BankHeader &
  get_obj_header()
  {
    return obj_header;
  }

  void
  set_obj_header(const BankHeader &objHeader)
  {
    obj_header = objHeader;
  }

  //! record bankwapper to obj_header and obj_classname
  void
  set_obj_info(const PgPostCalBank *bankwapper);

 private:
  //! header as written in database table
  BankHeader database_header;
  //! header as written in the bank wrapper object
  BankHeader obj_header;
  //! class name of the PgPostCalBank object as in the database
  std::string obj_classname;
  //! storage of the original calibration bank (vector of PdbCalChan)
  PdbCalBank *bank;

  ClassDefOverride(PgPostBankBackupStorage, 1)
};

#endif /* PDBCAL_PG_PGPOSTBANKBACKUPSTORAGE_H */
