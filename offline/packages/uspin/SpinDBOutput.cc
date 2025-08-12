#include "SpinDBOutput.h"
#include "SpinDBContent.h"
#include "SpinDBContentv1.h"

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>

#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <sstream>

//////////////////////////////////////////////////////////

void SpinDBOutput::Initialize()
{
  // default get from spin table from spin DB
  db_name = "spinDB";
  user_name = "phnxrc";
  table_name = "spin";

  spin_cont_store.clear();
  return;
}

//////////////////////////////////////////////////////////

void SpinDBOutput::SetDBName(const std::string &dbname)
{
  // Only DB to read from is "spinDB"
  if (dbname == "spinDB")
  {
    std::cout << std::format(" Database name is changed from {} to {}", db_name, dbname) << std::endl;
    db_name = dbname;
  }
  else
  {
    std::cout << std::format(" Your input database name, {}, is invalid", dbname) << std::endl;
    std::cout << std::format(" No change, try it as default, {}", db_name) << std::endl;
  }

  return;
}

//////////////////////////////////////////////////////////

void SpinDBOutput::SetTableName(const std::string &tname)
{
  // Only table is "spin"
  if (tname == "spin")
  {
    std::cout << std::format(" Table_name is changed from {} to {}", table_name, tname) << std::endl;
    table_name = tname;
  }
  else
  {
    std::cout << std::format(" Your input table name, {}, is invalid", tname) << std::endl;
    std::cout << std::format(" No change, try it as default, {}", db_name) << std::endl;
  }
  return;
}

//////////////////////////////////////////////////////////

odbc::Connection *SpinDBOutput::ConnectDB()
{
  odbc::Connection *con = nullptr;
  try
  {
    con = odbc::DriverManager::getConnection(db_name, user_name, "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    return (nullptr);
  }

  return (con);
}

/////////////////////////////////////////////////////////////

int SpinDBOutput::PrintDBColumn()
{
  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  int ncol = rs->getMetaData()->getColumnCount();
  for (int icol = 0; icol < ncol; icol++)
  {
    int col_type = rs->getMetaData()->getColumnType(icol + 1);
    std::string col_type_name = rs->getMetaData()->getColumnTypeName(icol + 1);
    std::string col_name = rs->getMetaData()->getColumnName(icol + 1);
    std::cout << std::format("{:2} : {:2} {:7s} {}", icol, col_type, col_type_name, col_name) << std::endl;
  }

  delete rs;
  delete stmt;
  delete con;
  return (1);
}

///////////////////////////////////////////////////////////
int SpinDBOutput::PrintDBRawContent(int runnum)
{
  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber=" << runnum << " and is_default=TRUE;";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error : Can't find data for run {} ", runnum) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  int ncol = rs->getMetaData()->getColumnCount();
  for (int icol = 0; icol < ncol; icol++)
  {
    std::string col_name = rs->getMetaData()->getColumnName(icol + 1);
    std::string cont = rs->getString(col_name);
    std::cout << std::format("{:2} : {} = {}", icol, col_name, cont) << std::endl;
  }

  delete rs;
  delete stmt;
  delete con;
  return (1);
}

int SpinDBOutput::PrintDBRawContent(int runnum, int qa_level)
{
  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber=" << runnum << " and qa_level=" << qa_level << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error : Can't find data for run {} (with qa_level {})", runnum, qa_level) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  int ncol = rs->getMetaData()->getColumnCount();
  for (int icol = 0; icol < ncol; icol++)
  {
    std::string col_name = rs->getMetaData()->getColumnName(icol + 1);
    std::string cont = rs->getString(col_name);
    std::cout << std::format("{:2} : {} = {}", icol, col_name, cont) << std::endl;
  }

  delete rs;
  delete stmt;
  delete con;
  return (1);
}

///////////////////////////////////////////////////////////
int SpinDBOutput::CheckRunRow(int runnum)
{
  if (CheckRunRowStore(runnum) == 1)
  {
    return (1);
  }

  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    std::cout << "No Connection to the Database!";
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber=" << runnum << " and is_default=TRUE;";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  delete rs;
  delete stmt;
  delete con;

  return (1);
}

int SpinDBOutput::CheckRunRow(int runnum, int qa_level)
{
  if (CheckRunRowStore(runnum) == 1)
  {
    return (1);
  }

  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    std::cout << "No Connection to the Database!";
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber=" << runnum << " and qa_level=" << qa_level << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  delete rs;
  delete stmt;
  delete con;

  return (1);
}

///////////////////////////////////////////////////////////

int SpinDBOutput::CheckRunRowStore(int runnum)
{
  std::map<int, std::unique_ptr<SpinDBContent>>::iterator it = spin_cont_store.find(runnum);
  if (it == spin_cont_store.end())
  {
    return 0;
  }
  return 1;
}

/////////////////////////////////////////////////////////////
int SpinDBOutput::StoreDBContent(int run1, int run2)
{
  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber>" << run1 - 1;
  cmd << " and runnumber<" << run2 + 1 << " and is_default=TRUE;";

  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  while (rs->next() != 0)
  {
    std::unique_ptr<SpinDBContent> spin_cont(new SpinDBContentv1());
    GetDBContent(*spin_cont, rs);

    int runnum = spin_cont->GetRunNumber();
    std::map<int, std::unique_ptr<SpinDBContent>>::iterator it = spin_cont_store.find(runnum);
    if (it != spin_cont_store.end())
    {
      it->second = std::move(spin_cont);
    }
    else
    {
      spin_cont_store.insert(std::make_pair(runnum, std::move(spin_cont)));
    }
  }

  delete rs;
  delete stmt;
  delete con;

  return (1);
}

int SpinDBOutput::StoreDBContent(int run1, int run2, int qa_level)
{
  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber>" << run1 - 1;
  cmd << " and runnumber<" << run2 + 1 << " and qa_level=" << qa_level << ";";

  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  while (rs->next() != 0)
  {
    std::unique_ptr<SpinDBContent> spin_cont(new SpinDBContentv1());
    GetDBContent(*spin_cont, rs);

    int runnum = spin_cont->GetRunNumber();
    std::map<int, std::unique_ptr<SpinDBContent>>::iterator it = spin_cont_store.find(runnum);
    if (it != spin_cont_store.end())
    {
      it->second = std::move(spin_cont);
    }
    else
    {
      spin_cont_store.insert(std::make_pair(runnum, std::move(spin_cont)));
    }
  }

  delete rs;
  delete stmt;
  delete con;

  return (1);
}

/////////////////////////////////////////////////////////////////

void SpinDBOutput::ClearDBContent()
{
  spin_cont_store.clear();
  return;
}

//////////////////////////////////////////////////////////////
int SpinDBOutput::GetDBContent(SpinDBContent *&spin_cont, int runnum)
{
  std::map<int, std::unique_ptr<SpinDBContent>>::iterator it = spin_cont_store.find(runnum);
  if (it != spin_cont_store.end())
  {
    spin_cont = it->second.get();
    return 1;
  }

  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber=" << runnum << " and is_default=TRUE;";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error : Can't find data for run {}", runnum) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  std::unique_ptr<SpinDBContent> spin_cont_ptr(new SpinDBContentv1());
  GetDBContent(*spin_cont_ptr, rs);
  spin_cont_store[runnum] = std::move(spin_cont_ptr);
  CopyDBContent(*spin_cont_store[runnum], *spin_cont);

  delete rs;
  delete stmt;
  delete con;

  return (1);
}

int SpinDBOutput::GetDBContent(SpinDBContent *&spin_cont, int runnum, int qa_level)
{
  std::map<int, std::unique_ptr<SpinDBContent>>::iterator it = spin_cont_store.find(runnum);
  if (it != spin_cont_store.end())
  {
    spin_cont = it->second.get();
    return 1;
  }

  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select * from " << table_name << " where runnumber=" << runnum << " with qa_level=" << qa_level << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error : Can't find data for run {} (with qa_level {})", runnum, qa_level) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  std::unique_ptr<SpinDBContent> spin_cont_ptr(new SpinDBContentv1());
  GetDBContent(*spin_cont_ptr, rs);
  spin_cont_store[runnum] = std::move(spin_cont_ptr);
  CopyDBContent(*spin_cont_store[runnum], *spin_cont);

  delete rs;
  delete stmt;
  delete con;

  return (1);
}

/////////////////////////////////////////////////////////////

int SpinDBOutput::GetDBContentStore(SpinDBContent *&spin_cont, int runnum)
{
  std::map<int, std::unique_ptr<SpinDBContent>>::iterator it = spin_cont_store.find(runnum);
  if (it != spin_cont_store.end())
  {
    CopyDBContent(*it->second, *spin_cont);
    return 1;
  }

  std::cout << std::format("Error : Can't find row for run {}.", runnum) << std::endl;

  return 0;
}

//////////////////////////////////////////////////////////////
int SpinDBOutput::CopyDBContent(SpinDBContent &spin_cont, SpinDBContent &spin_cont_copy)
{
  int ncross = SpinDBContent::GetNCrossing();
  spin_cont_copy.SetRunNumber(spin_cont.GetRunNumber());
  spin_cont_copy.SetQALevel(spin_cont.GetQALevel());
  spin_cont_copy.SetFillNumber(spin_cont.GetFillNumber());
  spin_cont_copy.SetBadRunFlag(spin_cont.GetBadRunFlag());
  spin_cont_copy.SetCrossingShift(spin_cont.GetCrossingShift());

  for (int i = 0; i < ncross; i++)
  {
    float b_pol;
    float b_pol_err;
    float b_pol_sys;
    float y_pol;
    float y_pol_err;
    float y_pol_sys;
    spin_cont.GetPolarizationBlue(i, b_pol, b_pol_err, b_pol_sys);
    spin_cont.GetPolarizationYellow(i, y_pol, y_pol_err, y_pol_sys);
    spin_cont_copy.SetPolarizationBlue(i, b_pol, b_pol_err, b_pol_sys);
    spin_cont_copy.SetPolarizationYellow(i, y_pol, y_pol_err, y_pol_sys);
    spin_cont_copy.SetSpinPatternBlue(i, spin_cont.GetSpinPatternBlue(i));
    spin_cont_copy.SetSpinPatternYellow(i, spin_cont.GetSpinPatternYellow(i));
    spin_cont_copy.SetScalerMbdVertexCut(i, spin_cont.GetScalerMbdVertexCut(i));
    spin_cont_copy.SetScalerMbdNoCut(i, spin_cont.GetScalerMbdNoCut(i));
    spin_cont_copy.SetScalerZdcNoCut(i, spin_cont.GetScalerZdcNoCut(i));
    spin_cont_copy.SetBadBunchFlag(i, spin_cont.GetBadBunchFlag(i));
  }

  float asym_bf;
  float asymerr_bf;
  float asym_bb;
  float asymerr_bb;
  float asym_yf;
  float asymerr_yf;
  float asym_yb;
  float asymerr_yb;
  float phase_bf;
  float phaseerr_bf;
  float phase_bb;
  float phaseerr_bb;
  float phase_yf;
  float phaseerr_yf;
  float phase_yb;
  float phaseerr_yb;
  spin_cont.GetAsymBlueForward(asym_bf, asymerr_bf);
  spin_cont.GetAsymBlueBackward(asym_bb, asymerr_bb);
  spin_cont.GetAsymYellowForward(asym_yf, asymerr_yf);
  spin_cont.GetAsymYellowBackward(asym_yb, asymerr_yb);
  spin_cont.GetPhaseBlueForward(phase_bf, phaseerr_bf);
  spin_cont.GetPhaseBlueBackward(phase_bb, phaseerr_bb);
  spin_cont.GetPhaseYellowForward(phase_yf, phaseerr_yf);
  spin_cont.GetPhaseYellowBackward(phase_yb, phaseerr_yb);

  spin_cont_copy.SetAsymBlueForward(asym_bf, asymerr_bf);
  spin_cont_copy.SetAsymBlueBackward(asym_bb, asymerr_bb);
  spin_cont_copy.SetAsymYellowForward(asym_yf, asymerr_yf);
  spin_cont_copy.SetAsymYellowBackward(asym_yb, asymerr_yb);
  spin_cont_copy.SetPhaseBlueForward(phase_bf, phaseerr_bf);
  spin_cont_copy.SetPhaseBlueBackward(phase_bb, phaseerr_bb);
  spin_cont_copy.SetPhaseYellowForward(phase_yf, phaseerr_yf);
  spin_cont_copy.SetPhaseYellowBackward(phase_yb, phaseerr_yb);

  spin_cont_copy.SetCrossAngle(spin_cont.GetCrossAngle());
  spin_cont_copy.SetCrossAngleStd(spin_cont.GetCrossAngleStd());
  spin_cont_copy.SetCrossAngleMin(spin_cont.GetCrossAngleMin());
  spin_cont_copy.SetCrossAngleMax(spin_cont.GetCrossAngleMax());

  return (1);
}

int SpinDBOutput::GetDBContent(SpinDBContent &spin_cont, odbc::ResultSet *rs)
{
  int ncross = SpinDBContent::GetNCrossing();

  spin_cont.SetRunNumber(rs->getInt("runnumber"));
  spin_cont.SetQALevel(rs->getInt("qa_level"));
  spin_cont.SetFillNumber(rs->getInt("fillnumber"));
  spin_cont.SetBadRunFlag(rs->getInt("badrunqa"));
  spin_cont.SetCrossingShift(rs->getInt("crossingshift"));

  // float bpol[ncross], bpolerr[ncross], bpolsys[ncross], ypol[ncross], ypolerr[ncross], ypolsys[ncross];
  // int bpat[ncross], ypat[ncross], bad_bunch[ncross];
  // long long mbd_vtxcut[ncross], mbd_nocut[ncross];
  // long long zdc_nocut[ncross];

  float *bpol = new float[ncross];
  float *bpolerr = new float[ncross];
  float *bpolsys = new float[ncross];
  float *ypol = new float[ncross];
  float *ypolerr = new float[ncross];
  float *ypolsys = new float[ncross];

  int *bpat = new int[ncross];
  int *ypat = new int[ncross];
  int *bad_bunch = new int[ncross];

  long long *mbd_vtxcut = new long long[ncross];
  long long *mbd_nocut = new long long[ncross];
  long long *zdc_nocut = new long long[ncross];

  GetArray(rs, "polarblue", bpol, ncross);
  GetArray(rs, "polarblueerror", bpolerr, ncross);
  if (table_name == "spin")
  {
    GetArray(rs, "polarblueerrorsys", bpolsys, ncross);
  }

  GetArray(rs, "polaryellow", ypol, ncross);
  GetArray(rs, "polaryellowerror", ypolerr, ncross);
  if (table_name == "spin")
  {
    GetArray(rs, "polaryellowerrorsys", ypolsys, ncross);
  }
  else
  {
    for (int i = 0; i < ncross; i++)
    {
      ypolsys[i] = std::numeric_limits<float>::quiet_NaN();
    }
  }
  GetArray(rs, "spinpatternblue", bpat, ncross);
  GetArray(rs, "spinpatternyellow", ypat, ncross);
  GetArray(rs, "mbdvtx", mbd_vtxcut, ncross);
  GetArray(rs, "mbdns", mbd_nocut, ncross);
  GetArray(rs, "zdcns", zdc_nocut, ncross);
  GetArray(rs, "badbunchqa", bad_bunch, ncross);

  for (int i = 0; i < ncross; i++)
  {
    spin_cont.SetPolarizationBlue(i, bpol[i], bpolerr[i], bpolsys[i]);
    spin_cont.SetPolarizationYellow(i, ypol[i], ypolerr[i], ypolsys[i]);
    spin_cont.SetSpinPatternBlue(i, bpat[i]);
    spin_cont.SetSpinPatternYellow(i, ypat[i]);
    spin_cont.SetScalerMbdVertexCut(i, mbd_vtxcut[i]);
    spin_cont.SetScalerMbdNoCut(i, mbd_nocut[i]);
    spin_cont.SetScalerZdcNoCut(i, zdc_nocut[i]);
    spin_cont.SetBadBunchFlag(i, bad_bunch[i]);
  }

  spin_cont.SetAsymBlueForward(rs->getFloat("asymbf"), rs->getFloat("asymerrbf"));
  spin_cont.SetAsymBlueBackward(rs->getFloat("asymbb"), rs->getFloat("asymerrbb"));
  spin_cont.SetAsymYellowForward(rs->getFloat("asymyf"), rs->getFloat("asymerryf"));
  spin_cont.SetAsymYellowBackward(rs->getFloat("asymyb"), rs->getFloat("asymerryb"));
  spin_cont.SetPhaseBlueForward(rs->getFloat("phasebf"), rs->getFloat("phaseerrbf"));
  spin_cont.SetPhaseBlueBackward(rs->getFloat("phasebb"), rs->getFloat("phaseerrbb"));
  spin_cont.SetPhaseYellowForward(rs->getFloat("phaseyf"), rs->getFloat("phaseerryf"));
  spin_cont.SetPhaseYellowBackward(rs->getFloat("phaseyb"), rs->getFloat("phaseerryb"));

  spin_cont.SetCrossAngle(rs->getFloat("crossingangle"));
  spin_cont.SetCrossAngleStd(rs->getFloat("crossanglestd"));
  spin_cont.SetCrossAngleMin(rs->getFloat("crossanglemin"));
  spin_cont.SetCrossAngleMax(rs->getFloat("crossanglemax"));

  delete[] bpol;
  delete[] bpolerr;
  delete[] bpolsys;
  delete[] ypol;
  delete[] ypolerr;
  delete[] ypolsys;

  delete[] bpat;
  delete[] ypat;
  delete[] bad_bunch;

  delete[] mbd_vtxcut;
  delete[] mbd_nocut;
  delete[] zdc_nocut;
  return (1);
}

///////////////////////////////////////////////////////////////

int SpinDBOutput::GetArray(odbc::ResultSet *rs, const std::string &name, std::vector<std::string> &value) const
{
  std::string cvalue;
  try
  {
    cvalue = rs->getString(name);
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    return (ERROR_VALUE);
  }

  int length = cvalue.size();
  if (length)
  {
    if (cvalue.compare(0, 1, "{") != 0 || cvalue.compare(length - 1, 1, "}") != 0)
    {
      if (verbosity > 0)
      {
        std::cout << std::format("Error : Is this array? ({}), length = {}", name, length) << std::endl;
      }
      return 0;
    }
  }
  else
  {
    if (verbosity > 0)
    {
      std::cout << std::format("Error : Is this array? ({}), length = {}", name, length) << std::endl;
    }
    return 0;
  }

  cvalue = cvalue.substr(1, cvalue.size() - 1) + ",";

  //  int nvalue = 0;
  std::vector<std::string> value1;
  while (true)
  {
    size_t pos = cvalue.find(std::string(","));
    if (pos == std::string::npos)
    {
      break;
    }

    value1.push_back(cvalue.substr(0, pos));
    cvalue.erase(0, pos + 1);
    //    nvalue++;
  }

  std::stringstream cerror;
  cerror << SpinDBContent::GetErrorValue();
  for (int i = 0; i < (int) value.size(); i++)
  {
    if (i < (int) value1.size())
    {
      value[i] = value1[i];
    }
    else
    {
      value[i] = cerror.str();
    }
  }

  if (value1.size() != value.size())
  {
    std::cout << "Error : Number of values are inconsistent with expected.";
    std::cout << " (" << name << " : " << value1.size() << " != " << value.size() << ")" << std::endl;
    return (0);
  }

  return (1);
}

//////////////////////////////////////////////////////////////

int SpinDBOutput::GetArray(odbc::ResultSet *rs, const std::string &name, float *value, int nvalue)
{
  std::vector<std::string> svalue(nvalue, "");
  int ret = GetArray(rs, name, svalue);
  for (int i = 0; i < nvalue; i++)
  {
    value[i] = atof(svalue[i].c_str());
  }
  return (ret);
}

/////////////////////////////////////////////////////////////////

int SpinDBOutput::GetArray(odbc::ResultSet *rs, const std::string &name, unsigned int *value, int nvalue)
{
  std::vector<std::string> svalue(nvalue, "");
  int ret = GetArray(rs, name, svalue);
  for (int i = 0; i < nvalue; i++)
  {
    value[i] = (unsigned int) atoi(svalue[i].c_str());
  }
  return (ret);
}

/////////////////////////////////////////////////////////////////

int SpinDBOutput::GetArray(odbc::ResultSet *rs, const std::string &name, int *value, int nvalue)
{
  std::vector<std::string> svalue(nvalue, "");
  int ret = GetArray(rs, name, svalue);
  for (int i = 0; i < nvalue; i++)
  {
    value[i] = atoi(svalue[i].c_str());
  }
  return (ret);
}

/////////////////////////////////////////////////////////////////

int SpinDBOutput::GetArray(odbc::ResultSet *rs, const std::string &name, long long *value, int nvalue)
{
  std::vector<std::string> svalue(nvalue, "");
  int ret = GetArray(rs, name, svalue);
  for (int i = 0; i < nvalue; i++)
  {
    value[i] = 0;
  }
  for (int i = 0; i < nvalue; i++)
  {
    std::istringstream iss(svalue[i]);
    iss >> value[i];
    // sscanf(svalue[i].c_str(),"%lld",&value[i]);
  }
  return (ret);
}

/////////////////////////////////////////////////////////////////

int SpinDBOutput::GetDefaultQA(int runnum)
{
  odbc::Connection *con = ConnectDB();
  if (con == nullptr)
  {
    return (0);
  }

  std::stringstream cmd;
  cmd << "select qa_level from " << table_name << " where runnumber=" << runnum << " and is_default=TRUE;";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << std::format("Error: {}.", e.getMessage()) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (ERROR_VALUE);
  }

  if (rs->next() == 0)
  {
    std::cout << std::format("Error : Can't find data for run {}", runnum) << std::endl;
    delete rs;
    delete stmt;
    delete con;
    return (0);
  }

  int default_qa_level = rs->getInt("qa_level");

  delete rs;
  delete stmt;
  delete con;

  return (default_qa_level);
}
