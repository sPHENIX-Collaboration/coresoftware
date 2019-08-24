/*!
 * \file PgPostBankBackupManager.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.11 $
 * \date $Date: 2014/11/24 17:46:49 $
 */

#include "PgPostBankBackupManager.h"

#include "PgPostApplication.h"
#include "PgPostBankBackupLog.h"
#include "PgPostBankBackupStorage.h"
#include "PgPostBankWrapper.h"
#include "PgPostCalBank.h"

#include <pdbcalbase/PdbCalChan.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/PdbClassMap.h>

#include <phool/PHTimeServer.h>
#include <phool/PHTimeStamp.h>
#include <phool/PHTimer.h>

#include <RDBC/TSQL.h>
#include <RDBC/TSQLConnection.h>
#include <RDBC/TSQLPreparedStatement.h>
#include <RDBC/TSQLResultSet.h>
#include <RDBC/TSQLStatement.h>

#include <TBuffer.h>
#include <TBufferFile.h>
#include <TCollection.h>                 // for TIter
#include <TDirectory.h>                  // for gDirectory, TDirectory (ptr ...
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>                     // for TObject, TObject::kWriteDelete
#include <TString.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

using namespace std;

PgPostBankBackupManager::PgPostBankBackupManager(const std::string &Tag)
  : verbosity(0)
  , tag(Tag)
{
  if (tag.length() == 0)
  {
    tag = "NO_TAG";
  }

  cout
      << "PgPostBankBackupManager::PgPostBankBackupManager - use backup log tag "
      << tag << endl;
}

PgPostBankBackupManager::~PgPostBankBackupManager()
{
}

void PgPostBankBackupManager::deleteSQLStatement(TSQLStatement *stmt)
{
  if (!stmt)
  {
    cout << "PgPostBankBackupManager::deleteSQLStatement - WARNING - "
         << " empty pointer" << endl;
    return;
  }

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::deleteSQLStatement - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }
  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::deleteSQLStatement - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  TList *sl = con->GetListOfStatements();
  if (!sl)
    return;

  TObject *obj = sl->Remove(stmt);
  if (!obj)
  {
    cout << "PgPostBankBackupManager::deleteSQLStatement - ERROR - "
         << " Cannot find statement in TSQLConnection" << endl;
    exit(1);
  }

  delete stmt;
  //  stmt = nullptr;
}

void PgPostBankBackupManager::deleteODBCPreparedStatement(
    TSQLPreparedStatement *stmt)
{
  if (!stmt)
  {
    cout
        << "PgPostBankBackupManager::deleteODBCPreparedStatement - WARNING - "
        << " empty pointer" << endl;
    return;
  }

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::deleteODBCPreparedStatement - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }
  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::deleteODBCPreparedStatement - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  TList *sl = con->GetListOfStatements();
  if (!sl)
    return;

  // inherance check
  TObject *o = dynamic_cast<TObject *>(stmt);
  assert(o);

  TObject *obj = sl->Remove(o);
  if (!obj)
  {
    cout << "PgPostBankBackupManager::deleteODBCPreparedStatement - ERROR - "
         << " Cannot find statement in TSQLConnection" << endl;
    exit(1);
  }
  assert(o == obj);

  //  cout << "PgPostBankBackupManager::deleteODBCPreparedStatement - delete the object" << endl;
  delete stmt;
  //  stmt = nullptr;
}

PgPostBankBackupStorage *
PgPostBankBackupManager::SQLResultSet2BackupStorage(TSQLResultSet *rs,
                                                    const std::string &table_name)
{
  if (!rs)
  {
    cout
        << "PgPostBankBackupManager::SQLResultSet2BackupStorage - Error - null TSQLResultSet"
        << endl;
    return nullptr;
  }

  if (rs->GetRow() == 0)
  {
    cout
        << "PgPostBankBackupManager::SQLResultSet2BackupStorage - Error - invalid TSQLResultSet with GetRow = 0"
        << endl;
    return nullptr;
  }

  const int rid = rs->GetInt(8);
  std::unique_ptr<PgPostCalBank> bw(
      dynamic_cast<PgPostCalBank *>(rs->GetObject(7)));
  assert(bw.get());
  int length = 0;
  string pdbcalchan_classname;
  string pgpostcalbank_classname;
  PdbCalBank *bank_orig = nullptr;

  PdbClassMap<PdbCalBank> *classMap = PdbClassMap<PdbCalBank>::instance();
  if (string(bw->ClassName()) == string("PgPostBankWrapper"))
  {
    if (verbosity >= 2)
      cout
          << "PgPostBankBackupManager::SQLResultSet2BackupStorage - Processing PgPostBankWrapper "
          << endl;
    bank_orig = (static_cast<PgPostBankWrapper *>(bw.get()))->getBank();
    assert(bank_orig);
  }
  else if (string(bw->ClassName()) == string("PgPostCalBank"))
  {
    cout << "PgPostBankBackupManager::SQLResultSet2BackupStorage - WARNING - "
         << "empty PgPostCalBank object in database, table " << table_name
         << " where rid = " << rid << endl;
    bank_orig = dynamic_cast<PdbCalBank *>(bw.get());
    assert(bank_orig);
  }
  else if (classMap->find(bw->ClassName()) != classMap->end())
  {
    cout << "PgPostBankBackupManager::SQLResultSet2BackupStorage - WARNING - "
         << "Direct stream of " << bw->ClassName()
         << " object without wrapper in database, table " << table_name
         << " where rid = " << rid << endl;
    bank_orig = dynamic_cast<PdbCalBank *>(bw.get());
    assert(bank_orig);
  }
  else
  {
    cout << "PgPostBankBackupManager::SQLResultSet2BackupStorage - "
         << " ERROR - unknown object in database record in class "
         << bw->ClassName() << endl;
    exit(1);
  }

  // check bank_orig
  int unwrap_cnt = 1;
  while (bank_orig)
  {
    if (string(bank_orig->ClassName()) == string("PgPostBankWrapper"))
    {
      cout
          << "PgPostBankBackupManager::SQLResultSet2BackupStorage - WARNING - "
          << "PgPostBankWrapper object nested inside PgPostBankWrapper layer "
          << unwrap_cnt << ". Discard secondary wrapper: table "
          << table_name << " where rid = " << rid << endl;
      bank_orig = static_cast<PgPostBankWrapper *>(bank_orig)->getBank();
      assert(bank_orig);
    }
    else
    {
      if (verbosity >= 2)
        cout << "PgPostBankBackupManager::SQLResultSet2BackupStorage - "
             << bank_orig->ClassName() << " passed wrapper check" << endl;
      break;  // Done
    }

    unwrap_cnt++;
    if (unwrap_cnt > 10)
    {
      cout
          << "PgPostBankBackupManager::SQLResultSet2BackupStorage - Fatal Error - "
          << "PgPostBankWrapper object nested inside PgPostBankWrapper for "
          << unwrap_cnt << " layers: table " << table_name
          << " where rid = " << rid << endl;
      exit(1);
    }
  }

  // Check again
  assert(bank_orig);
  if (classMap->find(bank_orig->ClassName()) == classMap->end())
  {
    cout
        << "PgPostBankBackupManager::SQLResultSet2BackupStorage - Fatal Error - "
        << "Calibration bank object of " << bank_orig->ClassName()
        << " is not supported: table " << table_name << " where rid = " << rid
        << endl;
    exit(10);
  }

  // Go on to make copies
  length = bank_orig->getLength();
  pgpostcalbank_classname = bank_orig->ClassName();
  pdbcalchan_classname = string("Pdb") + getBankBaseName(pgpostcalbank_classname);
  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::SQLResultSet2BackupStorage - Bank "
         << pgpostcalbank_classname << " contains " << length << " "
         << pdbcalchan_classname << endl;

  PdbCalBank *bank_copy = static_cast<PdbCalBank *>(bank_orig->Clone());
  PgPostBankBackupStorage *bs = new PgPostBankBackupStorage(bank_copy);
  assert(bs);

  if (verbosity >= 40)
  {
    cout
        << "PgPostBankBackupManager::SQLResultSet2BackupStorage - obj header 1 "
        << endl;
    bs->get_obj_header().Print();
  }

  bs->get_database_header().setBankID(rs->GetInt(1));
  //      bs->get_database_header().setBankID2(rs->GetInt(1));
  bs->get_database_header().setInsertTime(rs->GetLong(2));
  bs->get_database_header().setStartValTime(rs->GetLong(3));
  bs->get_database_header().setEndValTime(rs->GetLong(4));
  bs->get_database_header().setDescription(rs->GetString(5).Data());
  bs->get_database_header().setUserName(rs->GetString(6).Data());
  bs->get_database_header().setTableName(table_name);
  bs->get_database_header().setRId(rid);

  if (verbosity >= 40)
  {
    cout
        << "PgPostBankBackupManager::SQLResultSet2BackupStorage - obj header 2 "
        << endl;
    bs->get_obj_header().Print();
  }

  bs->set_obj_info(bw.get());
  bs->format_name_title();

  if (verbosity >= 40)
  {
    cout
        << "PgPostBankBackupManager::SQLResultSet2BackupStorage - obj header 3 "
        << endl;
    bs->get_obj_header().Print();
  }

  return bs;
}

//! Build PgPostBankBackupStorage from one database entry
//! using code should delete the returned object
PgPostBankBackupStorage *
PgPostBankBackupManager::fetchBank(const std::string &bankName, int rid)
{
  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::fetchBank - start on fetching "
         << bankName << " ID " << rid << endl;

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::fetchBank - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }

  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::fetchBank - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  TSQLStatement *stmt = con->CreateStatement();
  std::ostringstream tem;
  //  std::ostringstream t2;

  // cout << bankID.getInternalValue() << endl;
  tem
      << "select bankid,inserttime,startvaltime,endvaltime,description,username,calibrations,rid from "
      << bankName << " where rid = " << rid;

  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::fetchBank - database exe : " << tem.str()
         << endl;

  //  std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
  //  if (!rs.get())
  TSQLResultSet *rs(stmt->ExecuteQuery(tem.str().c_str()));
  if (!rs)
  {
    cout << "PgPostBankBackupManager::fetchBank - ERROR - "
         << " Cannot get TSQLResultSet from ExecuteQuery, exiting" << endl;
    exit(1);
  }

  if (rs->Next())
  {
    PgPostBankBackupStorage *bs = SQLResultSet2BackupStorage(rs, bankName);

    if (Verbosity() >= 2)
      cout << "PgPostBankBackupManager::fetchBank - clear TSQLResultSet"
           << endl;

    delete rs;
    rs = nullptr;

    if (Verbosity() >= 2)
      cout << "PgPostBankBackupManager::fetchBank - clear SQLStatement List"
           << endl;
    deleteSQLStatement(stmt);

    return bs;
  }
  else
  {
    cout << "PgPostBankBackupManager::fetchBank - ERROR -  NO Bank found : "
         << tem.str() << std::endl;

    if (Verbosity() >= 2)
      cout << "PgPostBankBackupManager::fetchBank - clear TSQLResultSet"
           << endl;
    if (rs)
    {
      if (Verbosity() >= 2)
        cout << "PgPostBankBackupManager::fetchBank - clear TSQLResultSet"
             << endl;

      delete rs;
      rs = nullptr;
    }

    return nullptr;
  }

  return nullptr;
}

//! Use PgPostBankBackupStorage to restore the original entry in database
//! \param[in] bs PgPostBankBackupStorage object. This function will not own this object
//! \return true if success
bool PgPostBankBackupManager::commit(PgPostBankBackupStorage *bs)
{
  if (!bs)
  {
    cout << "PgPostBankBackupManager::commit - Error - "
         << " invalid input pointer" << endl;
    return false;
  }
  if (!bs->isValid())
  {
    cout << "PgPostBankBackupManager::commit - Error - "
         << " invalid input object" << endl;
    return false;
  }

  //  std::unique_ptr<PgPostCalBank> bw(bs->createBank());
  //  assert(bw.get());

  try
  {
    PgPostApplication *ap = PgPostApplication::instance();
    assert(ap);
    TSQLConnection *con = ap->getConnection();
    assert(con);

    // Check
    if (isRIdExist(bs->get_database_header().getTableName(),
                   bs->get_database_header().getRId()))
    {
      cout << "PgPostBankBackupManager::commit - Error - "
           << " RID "
           << bs->get_database_header().getRId()
           << " already exist in table "
           << bs->get_database_header().getTableName() << ". STOP!" << endl;
      return false;
    }

    // Commit
    ostringstream sqlcmd;
    sqlcmd << "insert into " << bs->get_database_header().getTableName()
           << "(bankid,inserttime,startvaltime,endvaltime,description,username,calibrations,rid) values (?,?,?,?,?,?,?,?);";

    if (verbosity >= 2)
      cout << "PgPostBankBackupManager::fetchBank - database commit : "
           << sqlcmd.str() << endl;
    TSQLPreparedStatement *pstmt = con->PrepareStatement(
        sqlcmd.str().c_str());
    //  std::unique_ptr<TSQLPreparedStatement> pstmt(
    //      con->PrepareStatement(sqlcmd.str().c_str()));

    //      deleteODBCPreparedStatement(pstmt);

    pstmt->SetInt(1, bs->get_database_header().getBankID());
    pstmt->SetLong(2, (bs->get_database_header().getInsertTime()).getTics());
    pstmt->SetLong(3,
                   (bs->get_database_header().getStartValTime()).getTics());
    pstmt->SetLong(4, (bs->get_database_header().getEndValTime()).getTics());
    pstmt->SetString(5, (bs->get_database_header().getDescription()).c_str());
    pstmt->SetString(6, (bs->get_database_header().getUserName()).c_str());
    std::unique_ptr<PgPostCalBank> bw(bs->createBank());
    assert(bw.get());
    pstmt->SetObject(7, bw.get());
    pstmt->SetInt(8, bs->get_database_header().getRId());
    int res = 0;
    res = pstmt->ExecuteUpdate();
    try
    {
      con->Commit();
    }
    catch (TSQLException &e)
    {
      cout << "PgPostBankBackupManager::commit - Error - "
           << " Exception caught during connection->commit()" << endl;
      cout << e.GetMessage() << endl;
      //      delete pstmt;
      return false;
    }
    //  delete pstmt;
    if (res == 0)
    {
      cout << "PgPostBankBackupManager::commit - Error - "
           << "DATABASE: commit to "
           << bs->get_database_header().getTableName() << " failed. "
           << "Make sure you commit to the master database " << endl;
      exit(1);
    }
    assert(res == 1);

    if (Verbosity() >= 2)
      cout << "PgPostBankBackupManager::fetchBank - clear SQLStatement List"
           << endl;
    deleteODBCPreparedStatement(pstmt);

    if (Verbosity() >= 1)
      cout << "PgPostBankBackupManager::commit - Committed " << res << " row"
           << endl;
  }
  catch (std::exception &e)
  {
    cout << "PgPostBankBackupManager::commit - Error - "
         << " Exception caught " << endl;
    cout << e.what() << endl;
    //      delete pstmt;
    return false;
  }

  return true;
}

//! TFile -> Database for official operation use
int PgPostBankBackupManager::commitAllBankfromTFile(const std::string &input_file)
{
  if (Verbosity() >= 1)
    cout << "PgPostBankBackupManager::commitAllBankfromTFile - Loading TFile "
         << input_file << endl;
  TFile *f = new TFile(input_file.c_str());

  if (!f)
  {
    cout
        << "PgPostBankBackupManager::commitAllBankfromTFile - ERROR - can not open TFile "
        << input_file << endl;
    return 0;
  }
  if (!f->IsOpen())
  {
    cout
        << "PgPostBankBackupManager::commitAllBankfromTFile - ERROR - can not open TFile "
        << input_file << endl;
    return 0;
  }

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::commitAllBankfromTFile - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }

  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::commitAllBankfromTFile - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  int commit_cnt = 0;
  int skip_cnt = 0;
  const int file_cnt = f->GetListOfKeys()->GetSize();

  PHTimeServer::timer timer_loop(
      PHTimeServer::get()->insert_new("Record loop"));
  PHTimeServer::timer timer_db(PHTimeServer::get()->insert_new("Database"));
  PHTimeServer::timer timer_file(PHTimeServer::get()->insert_new("File write"));
  PHTimeServer::timer timer_log(
      PHTimeServer::get()->insert_new("Log database"));

  bool first_read = true;
  string table_name;
  rid_list_t existing_rids;
  TSQLPreparedStatement *pstmt = nullptr;
  shared_ptr<PgPostBankBackupLog> bklog(
      static_cast<PgPostBankBackupLog *>(nullptr));

  TIter next(f->GetListOfKeys());
  TObject *o = nullptr;
  while ((o = next()))
  {
    if (Verbosity() >= 1)
    {
      if (((commit_cnt + skip_cnt) % 10000 == 1) || (Verbosity() >= 2))
      {
        cout
            << "PgPostBankBackupManager::commitAllBankfromTFile - processing "
            << table_name << ": commit/skip/total = " << commit_cnt << "/"
            << skip_cnt << "/" << file_cnt << " to " << input_file
            << endl;

        timer_loop.get()->print_stat();
        timer_db.get()->print_stat();
        timer_file.get()->print_stat();
        timer_log.get()->print_stat();
      }
    }

    timer_loop.get()->restart();

    TKey *key = dynamic_cast<TKey *>(o);
    assert(key);

    if (Verbosity() >= 2)
    {
      cout
          << "PgPostBankBackupManager::commitAllBankfromTFile - Processing "
          << key->GetName() << " ----------------------------------------"
          << endl;
    }

    try
    {
      timer_file.get()->restart();

      std::unique_ptr<PgPostBankBackupStorage>
      //          PgPostBankBackupStorage *
      bs(dynamic_cast<PgPostBankBackupStorage *>(key->ReadObj()));
      timer_file.get()->stop();

      if (!bs.get())
      {
        cout
            << "PgPostBankBackupManager::commitAllBankfromTFile - ERROR - can not read "
            << key->GetName() << " from " << input_file << endl;
        continue;
      }

      if (first_read)
      {
        first_read = false;

        table_name = bs->get_database_header().getTableName();

        existing_rids = PgPostBankBackupManager::getListOfRId(table_name);

        if (Verbosity() >= 1)
          cout
              << "PgPostBankBackupManager::commitAllBankfromTFile - wrtting table "
              << table_name << " from TFile " << input_file
              << ". TFile record size = " << file_cnt
              << " database size = " << existing_rids.size() << endl;

        ostringstream sqlcmd;
        sqlcmd << "insert into " << table_name
               << "(bankid,inserttime,startvaltime,endvaltime,description,username,calibrations,rid) values (?,?,?,?,?,?,?,?);";

        if (verbosity >= 1)
          cout
              << "PgPostBankBackupManager::commitAllBankfromTFile - database commit statement : "
              << sqlcmd.str() << endl;
        pstmt = con->PrepareStatement(sqlcmd.str().c_str());

        bklog = make_shared<PgPostBankBackupLog>(table_name, tag);
        bklog->Init();
      }  // if (first_read)

      assert(table_name.length() > 0);
      assert(pstmt);
      assert(bklog.get());

      if (table_name != bs->get_database_header().getTableName())
      {
        cout
            << "PgPostBankBackupManager::commitAllBankfromTFile - ERROR - inconsistent table name for "
            << key->GetName() << " from " << input_file << ", expect "
            << table_name << endl;

        continue;
      }

      const int rid = bs->get_database_header().getRId();
      const bool existing_rid = std::find(existing_rids.begin(),
                                          existing_rids.end(), rid) != existing_rids.end();
      if (existing_rid)
      {
        if (Verbosity() >= 2)
        {
          cout
              << "PgPostBankBackupManager::commitAllBankfromTFile - existing records "
              << key->GetName() << " in database" << endl;
        }

        timer_log.get()->restart();
        bklog->Log(rid, PgPostBankBackupLog::kOptFile2Db_Skip);
        timer_log.get()->stop();
        skip_cnt++;

      }  // existing records
      else
      {  // new records

        if (Verbosity() >= 2)
        {
          cout
              << "PgPostBankBackupManager::commitAllBankfromTFile - submit new record "
              << key->GetName() << " to database" << endl;
        }

        timer_db.get()->restart();
        pstmt->SetInt(1, bs->get_database_header().getBankID());
        pstmt->SetLong(2,
                       (bs->get_database_header().getInsertTime()).getTics());
        pstmt->SetLong(3,
                       (bs->get_database_header().getStartValTime()).getTics());
        pstmt->SetLong(4,
                       (bs->get_database_header().getEndValTime()).getTics());
        pstmt->SetString(5,
                         (bs->get_database_header().getDescription()).c_str());
        pstmt->SetString(6,
                         (bs->get_database_header().getUserName()).c_str());
        std::unique_ptr<PgPostCalBank> bw(bs->createBank());
        assert(bw.get());
        pstmt->SetObject(7, bw.get());
        pstmt->SetInt(8, bs->get_database_header().getRId());
        int res = 0;
        res = pstmt->ExecuteUpdate();
        if (res != 1)
        {
          cout
              << "PgPostBankBackupManager::commitAllBankfromTFile - Error - "
              << "DATABASE: commit to "
              << bs->get_database_header().getTableName()
              << " failed with ExecuteUpdate()=" << res
              << ". Make sure you commit to a writable database "
              << endl;

          timer_log.get()->restart();
          bklog->Log(rid,
                     (PgPostBankBackupLog::enu_ops)(PgPostBankBackupLog::kOptFile2Db * PgPostBankBackupLog::kOptFailed));
          timer_log.get()->stop();
        }
        else
        {
          try
          {
            con->Commit();
          }
          catch (TSQLException &e)
          {
            cout
                << "PgPostBankBackupManager::commitAllBankfromTFile - Error - "
                << " Exception caught during connection->commit()"
                << endl;
            cout << e.GetMessage() << endl;
            //      delete pstmt;
            exit(1);
          }
          timer_log.get()->restart();
          bklog->Log(rid, PgPostBankBackupLog::kOptFile2Db);
          timer_log.get()->stop();

          commit_cnt++;
        }

        timer_db.get()->stop();
      }

      timer_loop.get()->stop();

    }  // try
    catch (std::exception &e)
    {
      cout << "PgPostBankBackupManager::commit - Error - "
           << " unhandled exception caught when process record "
           << key->GetName() << ":" << endl;
      cout << e.what() << endl;
      //      delete pstmt;
    }

  }  // main loop of records

  f->Close();

  if (Verbosity() >= 1)
  {
    if ((commit_cnt + skip_cnt) % 10000 == 1)
    {
      cout << "PgPostBankBackupManager::commitAllBankfromTFile - Done "
           << table_name << ": commit/skip/total = " << commit_cnt << "/"
           << skip_cnt << "/" << file_cnt << " to " << input_file << endl;

      timer_loop.get()->print_stat();
      timer_db.get()->print_stat();
      timer_file.get()->print_stat();
      timer_log.get()->print_stat();
    }
  }

  return commit_cnt;
}

//! Database -> TFile
int PgPostBankBackupManager::fetchBank2TFile(const std::string &bankName,
                                             const rid_list_t &rid_list, const std::string &out_file)
{
  TDirectory *gd = gDirectory;
  TFile f(out_file.c_str(), "update");

  if (f.IsZombie())
  {
    cout << "PgPostBankBackupManager::fetchBank2TFile - Error -"
         << " can not open file" << out_file;
    return 0;
  }

  int cnt = 0;
  for (rid_list_t::const_iterator it = rid_list.begin(); it != rid_list.end();
       ++it)
  {
    const int rid = *it;
    if (Verbosity() >= 1)
      cout << "PgPostBankBackupManager::fetchBank2TFile - Process "
           << bankName << " at record ID " << rid << " (" << cnt << "/"
           << rid_list.size() << ")" << endl;

    PgPostBankBackupStorage *bs = fetchBank(bankName, rid);

    if (!bs)
    {
      cout << "PgPostBankBackupManager::fetchBank2TFile - Error -"
           << " nullptr PgPostBankBackupStorage object produced" << endl;
    }
    else if (!bs->isValid())
    {
      cout << "PgPostBankBackupManager::fetchBank2TFile - Error -"
           << " invalid PgPostBankBackupStorage object produced" << endl;
      bs->Print();
    }
    else
    {
      bs->Write(bs->GetName(), TObject::kWriteDelete);

      delete bs;
      cnt++;
    }
  }

  f.Close();

  gDirectory = gd;

  cout << "PgPostBankBackupManager::fetchBank2TFile - " << cnt << " records in "
       << bankName << " saved to " << out_file << endl;
  return cnt;
}

//! Database -> TFile for all all records satisifying  record_selection
int PgPostBankBackupManager::fetchAllBank2TFile(const std::string &bankName,
                                                const std::string &record_selection, const std::string &out_file_base,
                                                int splitting_limit)
{
  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::fetchBank - start on fetching "
         << bankName << " for records criteria of [" << record_selection
         << "] to " << out_file_base << "*.root" << endl;

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }

  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  const int row_cnt = getTotalRowCount(bankName);
  if (row_cnt <= 0)
  {
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - WARNING - "
         << " no data in table " << bankName << endl;
    //      return row_cnt;
  }

  TSQLStatement *stmt = con->CreateStatement();
  assert(stmt);

  if (verbosity >= 1)
    cout
        << "PgPostBankBackupManager::fetchAllBank2TFile - reset output length to row count "
        << row_cnt << " for " << bankName << endl;
  stmt->SetMaxRows(row_cnt + 1);

  std::ostringstream tem;
  //  std::ostringstream t2;
  tem
      << "select bankid,inserttime,startvaltime,endvaltime,description,username,calibrations,rid from "
      << bankName;
  if (record_selection.length() > 0)
    tem << " where " << record_selection;
  tem << "   ORDER BY rid ASC";

  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - database exe : "
         << tem.str() << endl;

  TDirectory *gd = gDirectory;
  TFile f;
  int file_id = 0;
  int file_rec_cnt = 0;
  TString file_name;

  file_name.Form("%s_%04d.root", out_file_base.c_str(), file_id);
  f.Open(file_name, "recreate");
  if (f.IsZombie())
  {
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - Error -"
         << " can not open file" << file_name << endl;
    return 0;
  }
  if (verbosity >= 1)
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - writing new file "
         << file_name << endl;

  std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
  if (!rs || !rs.get())
  //  TSQLResultSet * rs(stmt->ExecuteQuery(tem.str().c_str()));
  //  if (!rs)
  {
    cout << "PgPostBankBackupManager::fetchAllBank2TFile - ERROR - "
         << " Cannot get TSQLResultSet from ExecuteQuery, exiting" << endl;
    exit(1);
  }

  //  if (verbosity >= 2)
  //    cout << "PgPostBankBackupManager::fetchAllBank2TFile - received "
  //        << rs->GetRowCount() << " rows"
  ////    <<", max row = "<<
  //        << endl;

  PgPostBankBackupLog bklog(bankName, tag);
  bklog.Init();

  int cnt = 0;

  PHTimeServer::timer timer_loop(
      PHTimeServer::get()->insert_new("Record loop"));
  PHTimeServer::timer timer_db(PHTimeServer::get()->insert_new("Database"));
  PHTimeServer::timer timer_file(PHTimeServer::get()->insert_new("File write"));
  PHTimeServer::timer timer_log(
      PHTimeServer::get()->insert_new("Log database"));

  while (rs->Next())
  {
    PgPostBankBackupLog::enu_ops success = PgPostBankBackupLog::kOptSuccess;

    timer_loop.get()->restart();

    if (verbosity >= 1)
    {
      if (cnt % 1000 == 0)
        cout << "PgPostBankBackupManager::fetchAllBank2TFile - processing "
             << bankName << ": " << cnt << "/" << row_cnt << " to "
             << file_name << endl;

      if (cnt % 10000 == 1)
      {
        timer_loop.get()->print_stat();
        timer_db.get()->print_stat();
        timer_file.get()->print_stat();
        timer_log.get()->print_stat();
      }
    }

    if (file_rec_cnt >= splitting_limit)
    {
      f.Close();
      ++file_id;
      file_rec_cnt = 0;
      file_name.Form("%s_%04d.root", out_file_base.c_str(), file_id);
      f.Open(file_name, "recreate");
      if (f.IsZombie())
      {
        cout << "PgPostBankBackupManager::fetchAllBank2TFile - Error -"
             << " can not open file" << file_name;
        return 0;
      }
      if (verbosity >= 1)
        cout
            << "PgPostBankBackupManager::fetchAllBank2TFile - writing new file  "
            << file_name << endl;
    }

    timer_db.get()->restart();
    const int rid = rs->GetInt(8);
    PgPostBankBackupStorage *bs = SQLResultSet2BackupStorage(rs.get(),
                                                             bankName);
    timer_db.get()->stop();

    if (bs)
    {
      timer_file.get()->restart();
      bs->Write(bs->GetName(), TObject::kWriteDelete);
      timer_file.get()->stop();

      delete bs;
      cnt++;
      file_rec_cnt++;
    }
    else
    {
      success = PgPostBankBackupLog::kOptFailed;

      cout << "PgPostBankBackupManager::fetchAllBank2TFile - Error - "
           << "invalid PgPostBankBackupStorage for row " << rs->GetRow()
           << endl;
    }

    timer_log.get()->restart();
    bklog.Log(rid,
              (PgPostBankBackupLog::enu_ops)(success * PgPostBankBackupLog::kOptBackup2File));
    timer_log.get()->stop();

    timer_loop.get()->stop();
  }

  f.Close();
  gDirectory = gd;

  cout << "PgPostBankBackupManager::fetchAllBank2TFile - " << cnt
       << " records in " << bankName << " saved to " << out_file_base << "*.root"
       << endl;
  return cnt;
}

std::string
PgPostBankBackupManager::getBankBaseName(const std::string &bank_classname)
{
  if (bank_classname.length() <= 10)
  {
    cout
        << "PgPostBankBackupManager::getBankBaseName - Error - unknown format for "
        << bank_classname << endl;
    return "invalid";
  }
  string prefix = bank_classname.substr(0, 6);
  string suffix = bank_classname.substr(bank_classname.length() - 4, 4);

  if (prefix != "PgPost")
  {
    cout
        << "PgPostBankBackupManager::getBankBaseName - Error - unknown prefix format for "
        << bank_classname << endl;
    return "invalid";
  }
  if (suffix != "Bank")
  {
    cout
        << "PgPostBankBackupManager::getBankBaseName - Error - unknown suffix format for "
        << bank_classname << " with " << suffix << endl;
    return "invalid";
  }

  return bank_classname.substr(6, bank_classname.length() - 6 - 4);
}

//! check whether an record with certain rid already there
bool PgPostBankBackupManager::isRIdExist(const std::string &bankName, int rid)
{
  PgPostApplication *ap = PgPostApplication::instance();
  assert(ap);
  TSQLConnection *con = ap->getConnection();
  assert(con);

  // Check
  std::ostringstream tem;
  //  std::ostringstream t2;

  tem << "select bankid, description,rid from " << bankName << " where rid = "
      << rid;

  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::isRIdExist - database check : "
         << tem.str() << endl;

  TSQLStatement *stmt = con->CreateStatement();
  std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
  if ((rs) && rs->Next())
  {
    return true;
  }

  return false;
}

//! total row counts
int PgPostBankBackupManager::getTotalRowCount(const std::string &bankName)
{
  PgPostApplication *ap = PgPostApplication::instance();
  assert(ap);
  TSQLConnection *con = ap->getConnection();
  assert(con);

  // Check
  std::ostringstream tem;
  //  std::ostringstream t2;

  tem << "select COUNT(*) AS NumberOfRows from " << bankName;

  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::getTotalRowCount - database get : "
         << tem.str() << endl;

  TSQLStatement *stmt = con->CreateStatement();
  std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
  if ((rs) && rs->Next())
  {
    return rs->GetInt(1);
  }

  return 0;
}

//! list rid in a bank
PgPostBankBackupManager::rid_list_t
PgPostBankBackupManager::getListOfRId(const string &bankName,
                                      const string &condition)
{
  rid_list_t l;

  PgPostApplication *ap = PgPostApplication::instance();
  assert(ap);
  TSQLConnection *con = ap->getConnection();
  assert(con);

  // Check
  std::ostringstream tem;
  //  std::ostringstream t2;

  tem << "select rid from " << bankName;

  if (condition.length() > 0)
    tem << " where  " << condition;

  if (verbosity >= 1)
    cout << "PgPostBankBackupManager::getListOfRId - database fetch : "
         << tem.str() << endl;

  TSQLStatement *stmt = con->CreateStatement();
  std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
  while ((rs) && rs->Next())
  {
    int rid = rs->GetInt(1);
    l.push_back(rid);
  }
  if (verbosity >= 1)
    cout << "PgPostBankBackupManager::getListOfRId - database fetch "
         << l.size() << " RIDs" << endl;

  return l;
}

//! dump databse to text
void PgPostBankBackupManager::dumpTable(const std::string &bankName,
                                        std::ostream &out)
{
  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::dumpTable - start on fetching "
         << bankName << " -> " << bankName << endl;

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::dumpTable - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }

  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::dumpTable - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  TSQLStatement *stmt = con->CreateStatement();
  std::ostringstream tem;
  //  std::ostringstream t2;

  // cout << bankID.getInternalValue() << endl;
  tem
      << "select bankid,inserttime,startvaltime,endvaltime,description,username,calibrations,rid from "
      << bankName << "  ORDER BY rid ASC";

  if (verbosity >= 2)
    cout << "PgPostBankBackupManager::dumpTable - database exe : " << tem.str()
         << endl;

  std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
  if (!rs.get())
  {
    cout << "PgPostBankBackupManager::dumpTable - ERROR - "
         << " Cannot get TSQLResultSet from ExecuteQuery, exiting" << endl;
    exit(1);
  }

  PgPostBankBackupLog bklog(bankName, tag);
  bklog.Init();
  int cnt = 0;
  while (rs->Next())
  {
    std::unique_ptr<PgPostCalBank> bw(
        dynamic_cast<PgPostCalBank *>(rs->GetObject(7)));
    assert(bw.get());

    cnt++;

    const int rid = rs->GetInt(8);
    out << "New Entry - " << bankName << endl;
    out << "\tRId = \t" << rid << endl;
    out << "\tBankID = \t" << (rs->GetInt(1)) << endl;
    out << "\tInsertTime = \t" << (rs->GetLong(2)) << endl;
    out << "\tStartValTime = \t" << (rs->GetLong(3)) << endl;
    out << "\tEndValTime = \t" << (rs->GetLong(4)) << endl;
    out << "\tDescription = \t" << (rs->GetString(5).Data()) << endl;
    out << "\tUserName = \t" << (rs->GetString(6).Data()) << endl;
    out << "\tN PdbCalChan = \t" << (bw->getLength()) << endl;

    ULong_t hash = 0x0;
    for (size_t i = 0; i < bw->getLength(); i++)
    {
      PdbCalChan &c = bw->getEntry(i);

      //          out << "\t\t" << c.ClassName() << "[" << i << "] = "
      //              << HashPdbCalChan(c) << endl;

      hash ^= HashPdbCalChan(c);

      if (verbosity >= 4)
      {
        cout << "PgPostBankBackupManager::dumpTable -  " << bankName
             << ":" << rid << ", " << c.GetUniqueID() << ", "
             << c.TestBits(0xFFFFFFFF) << " -> HASH " << HashPdbCalChan(c)
             << ", ";
        c.print();
      }
    }
    out << "\tHash[PdbCalChans] = \t" << hash << endl;

    bklog.Log(rid,
              (PgPostBankBackupLog::enu_ops)(PgPostBankBackupLog::kOptDump));
  }
  cout << "PgPostBankBackupManager::dumpTable - processed " << cnt << " records"
       << endl;
}

//! Hash the object
ULong_t
PgPostBankBackupManager::HashPdbCalChan(const PdbCalChan &c)
{
  const TObject *ptr = dynamic_cast<const TObject *>(&c);
  assert(ptr);
  std::unique_ptr<TBuffer> b(new TBufferFile(TBuffer::kWrite));

  b->WriteObject(ptr);

  //  cout <<"PgPostBankBackupManager::HashPdbCalChan - buffer dump with size "<<b->Length()<<endl;
  //  const char * arr = b->Buffer();
  //  for (int i = 0; i<b->Length(); i++ )
  //    {
  //      cout <<i<<"\t = 0x";
  //      cout <<std::hex<<static_cast<unsigned int>(arr[i])<<endl;
  //    }
  //  cout <<std::dec;

  return TString::Hash(b->Buffer(), b->Length());
}

int PgPostBankBackupManager::CleanTable(const std::string &bankName,
                                        const PHTimeStamp &min_save_time, const bool do_delete,
                                        const std::string &log_file_name)
{
  const bool do_log = log_file_name.size();

  if (verbosity >= 1)
  {
    cout << "PgPostBankBackupManager::CleanTable - start on fetching "
         << bankName << " -> " << bankName
         << " with preservation limit of T>" << min_save_time;
    if (do_log)
      cout << ". Log to file " << log_file_name;
    cout << endl;
  }

  fstream flog;
  if (do_log)
    flog.open(log_file_name.c_str(), ios_base::out);

  PgPostApplication *ap = PgPostApplication::instance();
  if (!ap)
  {
    cout << "PgPostBankBackupManager::CleanTable - ERROR - "
         << " PgPostApplication instance is nullptr, exiting" << endl;
    exit(1);
  }

  TSQLConnection *con = ap->getConnection();
  if (!con)
  {
    cout << "PgPostBankBackupManager::CleanTable - ERROR - "
         << " Cannot get TSQLConnection, exiting" << endl;
    exit(1);
  }

  vector<int> bankids;
  {  // get bank IDs
    TSQLStatement *stmt = con->CreateStatement();
    std::ostringstream tem;

    tem << "select bankid from " << bankName
        << " group by bankid order by bankid";

    if (verbosity >= 1)
      cout << "PgPostBankBackupManager::CleanTable - database exe : "
           << tem.str() << endl;

    std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
    if (!rs.get())
    {
      cout << "PgPostBankBackupManager::CleanTable - ERROR - "
           << " Cannot get TSQLResultSet from ExecuteQuery, exiting" << endl;
      exit(1);
    }

    while (rs->Next())
    {
      const int id = rs->GetInt(1);

      if (verbosity >= 2)
        cout << "\tBankID = \t" << id << endl;

      bankids.push_back(id);
    }

    if (verbosity >= 1)
      cout << "PgPostBankBackupManager::CleanTable -" << bankName
           << ": received " << bankids.size() << " bank IDs" << endl;
  }  // get bank IDs

  vector<int> rid_to_drop;
  int cnt_total = 0;
  for (vector<int>::const_iterator it = bankids.begin(); it != bankids.end();
       ++it)
  {
    const int bankid = *it;

    if (verbosity >= 1)
      cout << "PgPostBankBackupManager::CleanTable -" << bankName
           << ": process bankID " << bankid << ". Processed " << cnt_total
           << " records / delete " << rid_to_drop.size() << endl;

    TSQLStatement *stmt = con->CreateStatement();
    std::ostringstream tem;

    tem << "select inserttime,startvaltime,endvaltime,rid FROM " << bankName
        << " where bankid =  " << bankid << " and inserttime < "
        << min_save_time.getTics() << " order by inserttime desc, rid desc";

    if (verbosity >= 1)
      cout << "PgPostBankBackupManager::CleanTable - database exe : "
           << tem.str() << endl;

    std::unique_ptr<TSQLResultSet> rs(stmt->ExecuteQuery(tem.str().c_str()));
    if (!rs.get())
    {
      cout << "PgPostBankBackupManager::CleanTable - ERROR - "
           << " Cannot get TSQLResultSet from ExecuteQuery, exiting" << endl;
      exit(1);
    }

    typedef map<int, int> covered_period_t;
    covered_period_t covered_period;

    while (rs->Next())
    {
      cnt_total++;

      const int inserttime = rs->GetInt(1);
      int startvaltime = rs->GetInt(2);
      int endvaltime = rs->GetInt(3);
      const int rid = rs->GetInt(4);
      int drop_record = 0;

      if (verbosity >= 2)
        cout << "PgPostBankBackupManager::CleanTable -" << bankName
             << " : " << bankid << " : " << rid
             << " : process record with insert " << inserttime << " covers ("
             << startvaltime << ", " << endvaltime << ")" << endl;

      // move this part to DB query part
      //          if (inserttime >= min_save_time.getTics())
      //            {
      //              if (verbosity >= 2)
      //                {
      //                  cout << "PgPostBankBackupManager::CleanTable - " << bankName
      //                      << " : " << bankid << " : " << rid
      //                      << " : preserve record since insert time" << inserttime
      //                      << " > threshold " << min_save_time.getTics()
      //                      << ". Drop it \n";
      //                }
      //            }

      if (startvaltime > endvaltime)
      {
        if (verbosity >= 2)
        {
          cout << "PgPostBankBackupManager::CleanTable - ERROR -"
               << bankName << " : " << bankid << " : " << rid
               << ": wrong record start " << startvaltime << " > end "
               << endvaltime << ". Drop it \n";
          drop_record = 1;
        }
      }  // if (startvaltime > endvaltime)
      else
      {
        int last_starttime = -2;
        //              int last_endtime = -1;
        covered_period_t::iterator piter = covered_period.begin();
        while (piter != covered_period.end() && drop_record == 0)
        {
          int covered_startvaltime = piter->first;
          int covered_endvaltime = piter->second;

          //                  assert(last_starttime<last_endtime);
          assert(last_starttime < covered_startvaltime);

          const bool contain_startvaltime = startvaltime >= covered_startvaltime && startvaltime <= covered_endvaltime;
          const bool contain_endvaltime = endvaltime >= covered_startvaltime && endvaltime <= covered_endvaltime;

          if (contain_startvaltime && contain_endvaltime)
          {
            if (verbosity >= 2)
              cout << "PgPostBankBackupManager::CleanTable -"
                   << bankName << " : " << bankid << " : " << rid
                   << " : - fully covered record (" << startvaltime
                   << ", " << endvaltime << ") <-> ("
                   << covered_startvaltime << ", "
                   << covered_endvaltime << ")" << endl;

            drop_record = 2;
            ++piter;
          }  // if (contain_endvaltime)
          else if (

              contain_startvaltime && !contain_endvaltime)
          {
            if (verbosity >= 2)
              cout << "PgPostBankBackupManager::CleanTable -"
                   << bankName << " : " << bankid << " : " << rid
                   << " : - afterward merge record (" << startvaltime
                   << ", " << endvaltime << ") <-> ("
                   << covered_startvaltime << ", "
                   << covered_endvaltime << ")" << endl;

            assert(endvaltime > covered_endvaltime);

            //                    $covered_period{covered_startvaltime} = endvaltime;
            piter->second = endvaltime;
            startvaltime = covered_startvaltime;

            ++piter;

            //           $drop_record                           = 0;
          }  //else if (contain_startvaltime && !contain_endvaltime)
          else if (!contain_startvaltime && contain_endvaltime)
          {
            if (verbosity >= 2)

              cout << "PgPostBankBackupManager::CleanTable -"
                   << bankName << " : " << bankid << " : " << rid
                   << " : - forward merge record (" << startvaltime
                   << ", " << endvaltime << ") <-> ("
                   << covered_startvaltime << ", "
                   << covered_endvaltime << ")" << endl;

            assert(startvaltime < covered_startvaltime);

            //                    delete( covered_period{covered_startvaltime} );
            covered_period_t::iterator piter_tmp = piter;
            ++piter;
            covered_period.erase(piter_tmp);

            covered_startvaltime = startvaltime;
            covered_period[covered_startvaltime] = covered_endvaltime;
            endvaltime = covered_endvaltime;

            //           $drop_record                           = 0;

          }  // else if (!contain_startvaltime && (contain_endvaltime || endvaltime == covered_startvaltime - 1))
          else if (startvaltime < covered_startvaltime && endvaltime > covered_endvaltime)
          {
            if (verbosity >= 2)

              cout << "PgPostBankBackupManager::CleanTable -"
                   << bankName << " : " << bankid << " : " << rid
                   << " : - over merge record (" << startvaltime
                   << ", " << endvaltime << ") <-> ("
                   << covered_startvaltime << ", "
                   << covered_endvaltime << ")" << endl;

            //                    delete( $covered_period{covered_startvaltime} );
            covered_period_t::iterator piter_tmp = piter;
            ++piter;
            covered_period.erase(piter_tmp);

            covered_startvaltime = startvaltime;
            covered_endvaltime = endvaltime;
            covered_period[covered_startvaltime] = covered_endvaltime;

            //             $drop_record                           = 0;

          }  // if ( startvaltime < covered_startvaltime && endvaltime > covered_endvaltime )
          else
          {
            if (verbosity >= 2)

              cout << "PgPostBankBackupManager::CleanTable -"
                   << bankName << " : " << bankid << " : " << rid
                   << " : - not related records (" << startvaltime
                   << ", " << endvaltime << ") <-> ("
                   << covered_startvaltime << ", "
                   << covered_endvaltime << ")" << endl;

            ++piter;
          }  // else

          last_starttime = covered_startvaltime;
          //                  last_endtime = covered_endvaltime;
        }  // covered_period loop

      }  //if (startvaltime > endvaltime) - else

      if (drop_record > 0)
      {
        if (verbosity >= 2)
        {
          cout << "PgPostBankBackupManager::CleanTable -" << bankName
               << " : " << bankid << " : " << rid << " : - drop record"
               << endl;
        }

        rid_to_drop.push_back(rid);
      }
      else
      {
        if (verbosity >= 2)
        {
          cout << "PgPostBankBackupManager::CleanTable -" << bankName
               << " : " << bankid << " : " << rid << " : - keep record"
               << endl;
        }
        covered_period[startvaltime] = endvaltime;
      }  // if ($drop_record) - else

      if (verbosity >= 3)
      {
        cout << "PgPostBankBackupManager::CleanTable -" << bankName
             << " : " << bankid << " : " << rid
             << " : - print covered period";
        for (covered_period_t::iterator piter = covered_period.begin();
             piter != covered_period.end(); ++piter)
        {
          const int covered_startvaltime = piter->first;
          const int covered_endvaltime = piter->second;

          cout << ", (" << covered_startvaltime << ", "
               << covered_endvaltime << ")";
        }
        cout << endl;
      }  // if ( $verbosity >= 3 )

    }  // while (rs->Next()) - record loop

  }  // bankid loop

  // do the deleting
  {
    if (verbosity >= 1)
      cout << "PgPostBankBackupManager::CleanTable -" << bankName
           << ": will delete " << rid_to_drop.size() << " records." << endl;

    ostringstream sqlcmd;
    sqlcmd << "delete from " << bankName << " where rid = ?";

    if (verbosity >= 1)
      cout << "PgPostBankBackupManager::CleanTable - database commit : "
           << sqlcmd.str() << endl;

    TSQLPreparedStatement *pstmt = con->PrepareStatement(
        sqlcmd.str().c_str());

    PgPostBankBackupLog bklog(bankName, tag);
    bklog.Init();

    PHTimeServer::timer timer_loop(
        PHTimeServer::get()->insert_new("Record loop"));
    PHTimeServer::timer timer_db(PHTimeServer::get()->insert_new("Database"));
    PHTimeServer::timer timer_log(
        PHTimeServer::get()->insert_new("Log database"));

    int cnt = 0;
    for (vector<int>::const_iterator it = rid_to_drop.begin();
         it != rid_to_drop.end(); ++it)
    {
      timer_loop.get()->restart();

      const int rid = *it;

      if (verbosity >= 2)
        cout << "PgPostBankBackupManager::CleanTable -" << bankName
             << ": delete rid " << rid << ", execute flag (do_delete) = "
             << do_delete << endl;

      if (verbosity >= 1)
      {
        if (cnt % 100 == 0)
          cout << "PgPostBankBackupManager::CleanTable -  " << bankName
               << " delete : " << cnt << "/" << rid_to_drop.size() << endl;

        if (cnt % 100 == 1)
        {
          timer_loop.get()->print_stat();
          timer_db.get()->print_stat();
          timer_log.get()->print_stat();
        }
      }

      if (do_delete)
      {
        timer_db.get()->restart();
        pstmt->SetInt(1, rid);
        int res = 0;
        res = pstmt->ExecuteUpdate();
        try
        {
          con->Commit();
        }
        catch (TSQLException &e)
        {
          cout << "PgPostBankBackupManager::CleanTable - Error - "
               << " Exception caught during connection->commit()"
               << endl;
          cout << e.GetMessage() << endl;
          //      delete pstmt;
          exit(1);
        }
        timer_db.get()->stop();

        if (res == 0)
        {
          cout << "PgPostBankBackupManager::CleanTable - Error - "
               << "DATABASE: delete " << bankName << " failed. "
               << "Make sure you commit to the master database " << endl;
          exit(1);
        }
        assert(res == 1);

        if (do_log)
        {
          flog << "Deleted Entry - " << bankName << " RId " << rid
               << endl;
        }

        timer_log.get()->restart();
        bklog.Log(rid,
                  (PgPostBankBackupLog::enu_ops)(PgPostBankBackupLog::kOptDelete));
        timer_log.get()->stop();
      }

      cnt++;
      timer_loop.get()->stop();
    }

    timer_loop.get()->print_stat();
    timer_db.get()->print_stat();
    timer_log.get()->print_stat();
  }  // do the deleting

  return rid_to_drop.size();
}
