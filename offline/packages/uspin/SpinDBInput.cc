#include "SpinDBInput.h"

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/resultsetmetadata.h>

#include <boost/format.hpp>

#include <cstdlib>
#include <iostream>

const std::string SpinDBInput::DB_NAME = "spinDB_write";
const std::string SpinDBInput::TABLE_NAME = "spin";
const int SpinDBInput::ERROR_VALUE = -999;

SpinDBInput::~SpinDBInput()
{
  delete con;
}

///////////////////////////////////////////////////////

void SpinDBInput::Initialize()
{
  con = nullptr;
  run_check = -1;
  qa_check = -1;

  try
  {
    con = odbc::DriverManager::getConnection(DB_NAME, "phnxrc", "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    exit(1);
  }

  if (con != nullptr)
  {
    std::cout << (boost::format("Connected to %s DB.\n") % TABLE_NAME).str();
  }

  return;
}

////////////////////////////////////////////////////////

int SpinDBInput::IsConnected()
{
  if (con == nullptr)
  {
    std::cout << "Error : No connection to the DB.\n";
    return (0);
  }

  return (1);
}

////////////////////////////////////////////////////////////

int SpinDBInput::CheckRunRow(int runnum, int qa_level, const char *opt)
{
  if (std::string("simple") == opt && run_check == runnum && qa_check == qa_level)
  {
    return (1);
  }

  if (IsConnected() != 1)
  {
    return (ERROR_VALUE);
  }

  std::stringstream cmd;
  cmd << "select * from spin where runnumber=" << runnum << " and qa_level=" << qa_level << ";";
  odbc::Statement *stmt = con->createStatement();
  odbc::ResultSet *rs = nullptr;
  try
  {
    rs = stmt->executeQuery(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    delete rs;
    delete stmt;
    return (ERROR_VALUE);
  }

  int flag = 0;
  if (rs->next() != 0)
  {
    flag = 1;
    run_check = runnum;
    qa_check = qa_level;
  }

  delete rs;
  delete stmt;
  return (flag);
}

//////////////////////////////////////////////////////////

int SpinDBInput::CreateRunRow(int runnum, int qa_level)
{
  if (CheckRunRow(runnum, qa_level, "simple") == 1)
  {
    std::cout << (boost::format("SpinDBInput::CreateRunRow() Error : Row for run %1% with qa level %2% seems to exist, check again.\n") % runnum % qa_level).str();
    return (0);
  }

  std::stringstream cmd;
  cmd << "insert into spin (runnumber,qa_level) values(" << runnum << "," << qa_level << ");";
  odbc::Statement *stmt = con->createStatement();
  try
  {
    stmt->execute(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    return (ERROR_VALUE);
  }
  delete stmt;

  SpinDBContent spin_cont_temp;
  spin_cont_temp.SetRunNumber(runnum);
  spin_cont_temp.SetQALevel(qa_level);
  InitializeRunRow(spin_cont_temp);

  return (1);
}

///////////////////////////////////////////////////////////

int SpinDBInput::DeleteRunRow(int runnum, int qa_level)
{
  if (IsConnected() != 1)
  {
    return (0);
  }

  if (CheckRunRow(runnum, qa_level, "simple") != 1)
  {
    std::cout << (boost::format("SpinDBInput::DeleteRunRow() Error : Row for run %1% with qa level %2% seems not to exist, check again.\n") % runnum % qa_level).str();
    return (0);
  }

  std::stringstream cmd;
  cmd << "delete from spin where runnumber=" << runnum << " and qa_level=" << qa_level << ";";
  odbc::Statement *stmt = con->createStatement();
  try
  {
    stmt->execute(cmd.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    delete stmt;
    return (ERROR_VALUE);
  }

  run_check = -1;
  qa_check = -1;
  delete stmt;
  return (1);
}

///////////////////////////////////////////////////////////

int SpinDBInput::InitializeRunRow(SpinDBContent spin_cont)
{
  if (IsConnected() != 1)
  {
    return (0);
  }

  if (CheckRunRow(spin_cont.GetRunNumber(), spin_cont.GetQALevel(), "simple") != 1)
  {
    std::cout << (boost::format("SpinDBInput::UpdateDBContent() Error : Row for run %1% with qa level %2% seems not to exist, check again.\n") % spin_cont.GetRunNumber() % spin_cont.GetQALevel()).str();
    return (0);
  }

  int runnum = spin_cont.GetRunNumber();
  int qa_level = spin_cont.GetQALevel();

  int ncross = spin_cont.GetNCrossing();

  InitializeValue(runnum, qa_level, "fillnumber");
  InitializeValue(runnum, qa_level, "badrunqa");
  InitializeValue(runnum, qa_level, "crossingshift");

  InitializeArray(runnum, qa_level, "polarblue", ncross);
  InitializeArray(runnum, qa_level, "polarblueerror", ncross);
  InitializeArray(runnum, qa_level, "polarblueerrorsys", ncross);
  InitializeArray(runnum, qa_level, "polaryellow", ncross);
  InitializeArray(runnum, qa_level, "polaryellowerror", ncross);
  InitializeArray(runnum, qa_level, "polaryellowerrorsys", ncross);
  InitializeArray(runnum, qa_level, "spinpatternblue", ncross);
  InitializeArray(runnum, qa_level, "spinpatternyellow", ncross);
  InitializeArray(runnum, qa_level, "mbdvtx", ncross);
  InitializeArray(runnum, qa_level, "mbdns", ncross);
  InitializeArray(runnum, qa_level, "zdcns", ncross);
  InitializeArray(runnum, qa_level, "badbunchqa", ncross);

  InitializeValue(runnum, qa_level, "asymbf");
  InitializeValue(runnum, qa_level, "asymbb");
  InitializeValue(runnum, qa_level, "asymyf");
  InitializeValue(runnum, qa_level, "asymyb");
  InitializeValue(runnum, qa_level, "asymerrbf");
  InitializeValue(runnum, qa_level, "asymerrbb");
  InitializeValue(runnum, qa_level, "asymerryf");
  InitializeValue(runnum, qa_level, "asymerryb");

  InitializeValue(runnum, qa_level, "phasebf");
  InitializeValue(runnum, qa_level, "phasebb");
  InitializeValue(runnum, qa_level, "phaseyf");
  InitializeValue(runnum, qa_level, "phaseyb");
  InitializeValue(runnum, qa_level, "phaseerrbf");
  InitializeValue(runnum, qa_level, "phaseerrbb");
  InitializeValue(runnum, qa_level, "phaseerryf");
  InitializeValue(runnum, qa_level, "phaseerryb");

  InitializeValue(runnum, qa_level, "crossingangle");
  InitializeValue(runnum, qa_level, "crossanglestd");
  InitializeValue(runnum, qa_level, "crossanglemin");
  InitializeValue(runnum, qa_level, "crossanglemax");

  return (1);
}

///////////////////////////////////////////////////////////

int SpinDBInput::UpdateDBContent(SpinDBContent spin_cont)
{
  if (IsConnected() != 1)
  {
    return (0);
  }

  if (CheckRunRow(spin_cont.GetRunNumber(), spin_cont.GetQALevel(), "simple") != 1)
  {
    std::cout << (boost::format("SpinDBInput::UpdateDBContent() Error : Row for run %1% with qa level %2% seems not to exist, check again.\n") % spin_cont.GetRunNumber() % spin_cont.GetQALevel()).str();
    return (0);
  }

  int qa_level = spin_cont.GetQALevel();

  if (qa_level == ERROR_VALUE)
  {
    std::cout << "You did not set a qa_level.  Please do so with SpinDBContent::SetQALevel(int qa_level).  Check that the qa level you set does not exist for this run before trying again.\n";
    return (0);
  }

  int runnum = spin_cont.GetRunNumber();
  int ncross = spin_cont.GetNCrossing();
  int fillnum = spin_cont.GetFillNumber();
  int badrun = spin_cont.GetBadRunFlag();
  int xingshift = spin_cont.GetCrossingShift();
  if (fillnum != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "fillnumber", fillnum);
  }
  if (badrun != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "badrunqa", badrun);
  }
  if (xingshift != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "crossingshift", xingshift);
  }

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
  
  for (int i = 0; i < ncross; i++)
  {
    spin_cont.GetPolarizationBlue(i, bpol[i], bpolerr[i], bpolsys[i]);
    spin_cont.GetPolarizationYellow(i, ypol[i], ypolerr[i], ypolsys[i]);
    bpat[i] = spin_cont.GetSpinPatternBlue(i);
    ypat[i] = spin_cont.GetSpinPatternYellow(i);
    mbd_vtxcut[i] = spin_cont.GetScalerMbdVertexCut(i);
    mbd_nocut[i] = spin_cont.GetScalerMbdNoCut(i);
    zdc_nocut[i] = spin_cont.GetScalerZdcNoCut(i);
    bad_bunch[i] = spin_cont.GetBadBunchFlag(i);
  }

  bool cbpol = false;
  bool cbpolerr = false;
  bool cbpolsys = false;
  bool cypol = false;
  bool cypolerr = false;
  bool cypolsys = false;
  bool cbpat = false;
  bool cypat = false;
  bool cmbd_vtxcut = false;
  bool cmbd_nocut = false;
  bool czdc_nocut = false;
  bool cbad_bunch = false;

  for (int i = 0; i < ncross; i++)
  {
    if (bpol[i] != ERROR_VALUE)
    {
      cbpol = true;
    }
    if (bpolerr[i] != ERROR_VALUE)
    {
      cbpolerr = true;
    }
    if (bpolsys[i] != ERROR_VALUE)
    {
      cbpolsys = true;
    }
    if (ypol[i] != ERROR_VALUE)
    {
      cypol = true;
    }
    if (ypolerr[i] != ERROR_VALUE)
    {
      cypolerr = true;
    }
    if (ypolsys[i] != ERROR_VALUE)
    {
      cypolsys = true;
    }
    if (bpat[i] != ERROR_VALUE)
    {
      cbpat = true;
    }
    if (ypat[i] != ERROR_VALUE)
    {
      cypat = true;
    }
    if (mbd_vtxcut[i] != ERROR_VALUE)
    {
      cmbd_vtxcut = true;
    }
    if (mbd_nocut[i] != ERROR_VALUE)
    {
      cmbd_nocut = true;
    }
    if (zdc_nocut[i] != ERROR_VALUE)
    {
      czdc_nocut = true;
    }
    if (bad_bunch[i] != ERROR_VALUE)
    {
      cbad_bunch = true;
    }
  }

  if (cbpol)
  {
    UpdateArray(runnum, qa_level, "polarblue", bpol, ncross);
  }
  if (cbpolerr)
  {
    UpdateArray(runnum, qa_level, "polarblueerror", bpolerr, ncross);
  }
  if (cbpolsys)
  {
    UpdateArray(runnum, qa_level, "polarblueerrorsys", bpolsys, ncross);
  }
  if (cypol)
  {
    UpdateArray(runnum, qa_level, "polaryellow", ypol, ncross);
  }
  if (cypolerr)
  {
    UpdateArray(runnum, qa_level, "polaryellowerror", ypolerr, ncross);
  }
  if (cypolsys)
  {
    UpdateArray(runnum, qa_level, "polaryellowerrorsys", ypolsys, ncross);
  }
  if (cbpat)
  {
    UpdateArray(runnum, qa_level, "spinpatternblue", bpat, ncross);
  }
  if (cypat)
  {
    UpdateArray(runnum, qa_level, "spinpatternyellow", ypat, ncross);
  }
  if (cmbd_vtxcut)
  {
    UpdateArray(runnum, qa_level, "mbdvtx", mbd_vtxcut, ncross);
  }
  if (cmbd_nocut)
  {
    UpdateArray(runnum, qa_level, "mbdns", mbd_nocut, ncross);
  }
  if (czdc_nocut)
  {
    UpdateArray(runnum, qa_level, "zdcns", zdc_nocut, ncross);
  }
  if (cbad_bunch)
  {
    UpdateArray(runnum, qa_level, "badbunchqa", bad_bunch, ncross);
  }

  float a_bf, a_bb, a_yf, a_yb;
  float a_bf_err, a_bb_err, a_yf_err, a_yb_err;
  float p_bf, p_bb, p_yf, p_yb;
  float p_bf_err, p_bb_err, p_yf_err, p_yb_err;
  spin_cont.GetAsymBlueForward(a_bf, a_bf_err);
  spin_cont.GetAsymBlueBackward(a_bb, a_bb_err);
  spin_cont.GetAsymYellowForward(a_yf, a_yf_err);
  spin_cont.GetAsymYellowBackward(a_yb, a_yb_err);
  spin_cont.GetPhaseBlueForward(p_bf, p_bf_err);
  spin_cont.GetPhaseBlueBackward(p_bb, p_bb_err);
  spin_cont.GetPhaseYellowForward(p_yf, p_yf_err);
  spin_cont.GetPhaseYellowBackward(p_yb, p_yb_err);

  if (a_bf != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymbf", a_bf);
  }
  if (a_bf_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymerrbf", a_bf_err);
  }
  if (a_bb != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymbb", a_bb);
  }
  if (a_bb_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymerrbb", a_bb_err);
  }
  if (a_yf != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymyf", a_yf);
  }
  if (a_yf_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymerryf", a_yf_err);
  }
  if (a_yb != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymyb", a_yb);
  }
  if (a_yb_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "asymerryb", a_yb_err);
  }

  
  if (p_bf != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phasebf", p_bf);
  }
  if (p_bf_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phaseerrbf", p_bf_err);
  }
  if (p_bb != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phasebb", p_bb);
  }
  if (p_bb_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phaseerrbb", p_bb_err);
  }
  if (p_yf != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phaseyf", p_yf);
  }
  if (p_yf_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phaseerryf", p_yf_err);
  }
  if (p_yb != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phaseyb", p_yb);
  }
  if (p_yb_err != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "phaseerryb", p_yb_err);
  }



  float cross_angle = spin_cont.GetCrossAngle();
  float cross_angle_std = spin_cont.GetCrossAngleStd();
  float cross_angle_min = spin_cont.GetCrossAngleMin();
  float cross_angle_max = spin_cont.GetCrossAngleMax();

  if (cross_angle != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "crossingangle", cross_angle);
  }
  if (cross_angle_std != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "crossanglestd", cross_angle_std);
  }
  if (cross_angle_min != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "crossanglemin", cross_angle_min);
  }
  if (cross_angle_max != ERROR_VALUE)
  {
    UpdateValue(runnum, qa_level, "crossanglemax", cross_angle_max);
  }
  delete [] bpol;
  delete [] bpolerr;
  delete [] bpolsys;
  delete [] ypol;
  delete [] ypolerr;
  delete [] ypolsys;

  delete [] bpat;
  delete [] ypat;
  delete [] bad_bunch;

  delete [] mbd_vtxcut;
  delete [] mbd_nocut;
  delete [] zdc_nocut;

  return (1);
}

/////////////////////////////////////////////////////////

int SpinDBInput::UpdateValue(int runnum, int qa_level, const char *cmd)
{
  if (IsConnected() != 1)
  {
    return (0);
  }

  if (CheckRunRow(runnum, qa_level, "simple") != 1)
  {
    std::cout << (boost::format("SpinDBInput::UpdateDBContent() Error : Row for run %1% with qa level %2% seems not to exist, check again.\n") % runnum % qa_level).str();
    return (0);
  }

  odbc::Statement *stmt = con->createStatement();
  try
  {
    stmt->execute(cmd);
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    delete stmt;
    return (ERROR_VALUE);
  }
  delete stmt;
  std::cout << cmd << '\n';

  return (1);
}

////////////////////////////////////////////////////////

template <class T>
int SpinDBInput::UpdateValueTemp(int runnum, int qa_level, const char *name, T value)
{
  std::stringstream cmd;
  cmd << "update " << TABLE_NAME << " set " << name << "=" << value;
  cmd << " where runnumber=" << runnum << " and qa_level=" << qa_level << ";";
  return (UpdateValue(runnum, qa_level, cmd.str().c_str()));
}

//////////////////////////////////////////////////////////

int SpinDBInput::UpdateValue(int runnum, int qa_level, const char *name, int value)
{
  return (UpdateValueTemp(runnum, qa_level, name, value));
}

//////////////////////////////////////////////////////////

int SpinDBInput::UpdateValue(int runnum, int qa_level, const char *name, float value)
{
  return (UpdateValueTemp(runnum, qa_level, name, value));
}

//////////////////////////////////////////////////////////

int SpinDBInput::UpdateValue(int runnum, int qa_level, const char *name, double value)
{
  return (UpdateValueTemp(runnum, qa_level, name, value));
}

//////////////////////////////////////////////////////////

template <class T>
int SpinDBInput::UpdateArrayTemp(int runnum, int qa_level, const char *name, T *value, int nvalue)
{
  std::stringstream cmd;
  cmd << "update " << TABLE_NAME << " set " << name << "='{";
  for (int i = 0; i < nvalue; i++)
  {
    cmd << value[i];
    if (i < nvalue - 1)
    {
      cmd << ",";
    }
  }
  cmd << "}' where runnumber=" << runnum << " and qa_level=" << qa_level << ";";
  return (UpdateValue(runnum, qa_level, cmd.str().c_str()));
}

/////////////////////////////////////////////////////////////////

int SpinDBInput::UpdateArray(int runnum, int qa_level, const char *name, int *value, int nvalue)
{
  return (UpdateArrayTemp(runnum, qa_level, name, value, nvalue));
}

/////////////////////////////////////////////////////////

int SpinDBInput::UpdateArray(int runnum, int qa_level, const char *name, float *value, int nvalue)
{
  return (UpdateArrayTemp(runnum, qa_level, name, value, nvalue));
}

/////////////////////////////////////////////////////////

int SpinDBInput::UpdateArray(int runnum, int qa_level, const char *name, double *value, int nvalue)
{
  return (UpdateArrayTemp(runnum, qa_level, name, value, nvalue));
}

/////////////////////////////////////////////////////////

int SpinDBInput::UpdateArray(int runnum, int qa_level, const char *name, unsigned int *value, int nvalue)
{
  return (UpdateArrayTemp(runnum, qa_level, name, value, nvalue));
}

/////////////////////////////////////////////////////////

int SpinDBInput::UpdateArray(int runnum, int qa_level, const char *name, long long *value, int nvalue)
{
  return (UpdateArrayTemp(runnum, qa_level, name, value, nvalue));
}

//////////////////////////////////////////////////////////

int SpinDBInput::InitializeValue(int runnum, int qa_level, const char *name)
{
  return (UpdateValue(runnum, qa_level, name, ERROR_VALUE));
}

//////////////////////////////////////////////////////////

int SpinDBInput::InitializeArray(int runnum, int qa_level, const char *name, int nvalue)
{
  int *ERROR_ARRAY = new int[nvalue];
  for (int i = 0; i < nvalue; i++)
  {
    ERROR_ARRAY[i] = ERROR_VALUE;
  }
  int iret = UpdateArray(runnum, qa_level, name, ERROR_ARRAY, nvalue);
  delete [] ERROR_ARRAY;
  return (iret);
}

////////////////////////////////////////////////////////////
int SpinDBInput::SetDefaultQA(SpinDBContent spin_cont)
{
  if (IsConnected() != 1)
  {
    return (0);
  }


  int qa_level = spin_cont.GetQALevel();

  if (qa_level == ERROR_VALUE)
  {
    std::cout << "You did not set a qa_level.  Please do so with SpinDBContent::SetQALevel(int qa_level).  Check that the qa level you set does not exist for this run before trying again.\n";
    return (0);
  }

  int runnum = spin_cont.GetRunNumber();

  if (CheckRunRow(runnum, qa_level, "simple") != 1)
  {
    std::cout << (boost::format("SpinDBInput::DeleteRunRow() Error : Row for run %1% with qa level %2% seems not to exist, check again.\n") % runnum % qa_level).str();
    return (0);
  }

  std::stringstream cmd1;
  cmd1 << " update spin set is_default = FALSE where runnumber=" << runnum << " and is_default = TRUE;";
  odbc::Statement *stmt1 = con->createStatement();
  try
  {
    stmt1->execute(cmd1.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    delete stmt1;
    return (ERROR_VALUE);
  }


  std::stringstream cmd2;
  cmd2 << " update spin set is_default = TRUE where runnumber=" << runnum << " and qa_level = " << qa_level << ";";
  odbc::Statement *stmt2 = con->createStatement();
  try
  {
    stmt2->execute(cmd2.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << (boost::format("Error: %s.\n") % e.getMessage().c_str()).str();
    delete stmt2;
    return (ERROR_VALUE);
  }

  run_check = -1;
  qa_check = -1;
  delete stmt1;
  delete stmt2;
  return (1);
}
