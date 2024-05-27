#include "XingShiftCal.h"

#include <cdbobjects/CDBTTree.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TCanvas.h>
#include <TH1.h>

#include <boost/format.hpp>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>
#include <odbc++/types.h>

#include <iostream>

XingShiftCal::XingShiftCal(const std::string &name, const int poverwriteSpinEntry)
  : SubsysReco(name)
  , overwriteSpinEntry(poverwriteSpinEntry)
{
  // overwriteSpinEntry = poverwriteSpinEntry;
  nevt = 0;
  
  for (auto &scalercount : scalercounts)
  {
    for (unsigned long &j : scalercount)
    {
      j = 0;
    }
  }
  std::cout << "XingShiftCal::XingShiftCal(const std::string &name) Calling ctor" << std::endl;
}

XingShiftCal::~XingShiftCal()
{
  std::cout << "XingShiftCal::~XingShiftCal() Calling dtor" << std::endl;
}

int XingShiftCal::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "XingShiftCal::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int XingShiftCal::InitRun(PHCompositeNode * /*topNode*/)
{
  // std::cout << "XingShiftCal::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;

  recoConsts *rc = recoConsts::instance();
  runnumber = rc->get_IntFlag("RUNNUMBER");
  return Fun4AllReturnCodes::EVENT_OK;
}

int XingShiftCal::process_event(PHCompositeNode *topNode)
{
  if (done == 1)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  // std::cout << "XingShiftCal::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");

  if (evt->getEvtType() == BEGRUNEVENT)
  {
    //================ BeginRunEvent packets ================//
    pBluePol = evt->getPacket(packet_BLUEPOL);
    pYellPol = evt->getPacket(packet_YELLPOL);

    pBlueIntPattern = evt->getPacket(packet_BLUEINTPATTERN);
    pYellIntPattern = evt->getPacket(packet_YELLINTPATTERN);
    pBluePolPattern = evt->getPacket(packet_BLUEPOLPATTERN);
    pYellPolPattern = evt->getPacket(packet_YELLPOLPATTERN);

    pBlueFillNumber = evt->getPacket(packet_BLUEFILLNUMBER);
    pYellFillNumber = evt->getPacket(packet_YELLFILLNUMBER);
    //=======================================================//

    //========= Get beam polarizations ==============//
    polBlue = -999;
    polBlueErr = -999;
    if (pBluePol)
    {
      polBlue = pBluePol->iValue(0) / 10000.0;
      polBlueErr = pBluePol->iValue(1) / 10000.0;
      delete pBluePol;
    }

    polYellow = -999;
    polYellowErr = -999;
    if (pYellPol)
    {
      polYellow = pYellPol->iValue(0) / 10000.0;
      polYellowErr = pYellPol->iValue(1) / 10000.0;
      delete pYellPol;
    }
    //==========================================================//

    //============== Get intended spin patterns from buckets ==============//
    // there are 360 buckets for 120 bunches

    if (pBlueIntPattern && pBluePolPattern)
    {
      for (int i = 0; i < 360; i += 3)
      {
	blueFillPattern[i / 3] = pBlueIntPattern->iValue(i);
        if (pBlueIntPattern->iValue(i))
        {
          blueSpinPattern[i / 3] = pBluePolPattern->iValue(i);
        }
        else
        {
          blueSpinPattern[i / 3] = 10;
        }
      }
      delete pBlueIntPattern;
      delete pBluePolPattern;
    }

    if (pYellIntPattern && pYellPolPattern)
    {
      for (int i = 0; i < 360; i += 3)
      {
	yellFillPattern[i / 3] = pYellIntPattern->iValue(i);
        if (pYellIntPattern->iValue(i))
        {
          yellSpinPattern[i / 3] = pYellPolPattern->iValue(i);
        }
        else
        {
          yellSpinPattern[i / 3] = 10;
        }
      }
      delete pYellIntPattern;
      delete pYellPolPattern;
    }
    //=======================================================================//

    //============== Set fill number histogram ==============//
    fillnumberBlue = 0;
    fillnumberYellow = 0;
    if (pBlueFillNumber)
    {
      fillnumberBlue = pBlueFillNumber->iValue(0);
      delete pBlueFillNumber;
    }
    if (pYellFillNumber)
    {
      fillnumberYellow = pYellFillNumber->iValue(0);
      delete pYellFillNumber;
    }
    //=======================================================//
  }
  else if (evt->getEvtType() == DATAEVENT)
  {
    p = evt->getPacket(packet_GL1);
    int bunchnr = p->lValue(0, "BunchNumber");
    for (int i = 0; i < NTRIG; i++)
    {
      // 2nd arg of lValue: 0 is raw trigger count, 1 is live trigger count, 2 is scaled trigger count
      long gl1pscaler = p->lValue(i, "GL1PLIVE");
      scalercounts[i][bunchnr] = gl1pscaler;
    }
    delete p;
  }

  if (nevt > threshold)
  {
    Calibrate();
    if (success)
    {
      done = 1;
    }
    else
    {
      threshold += threshold;
    }
  }

  if (nevt > evtcap)
  {
    done = 1;
  }

  nevt++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int XingShiftCal::ResetEvent(PHCompositeNode * /*topNode*/)
{
  // std::cout << "XingShiftCal::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
/*
int XingShiftCal::EndRun(const int runnumber)
{
  std::cout << "XingShiftCal::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
*/
int XingShiftCal::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "XingShiftCal::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  if (done)
  {
    // if (!success){
    Calibrate(1);
    //}
    // CommitConstantsToFile();
  }
  else
  {
    success = false;
    std::cout << "Not enough statistics. Did not calibrate." << std::endl;
  }

  const std::string cdbfname = (boost::format("SPIN-%08d_crossingshiftCDBTTree.root") % runnumber).str();
  WriteToCDB(cdbfname);
  CommitToSpinDB();

  if (commitSuccessCDB)
  {
    std::cout << "Commit to CDB : SUCCESS" << std::endl;
  }
  else
  {
    std::cout << "Commit to CDB : FAILURE" << std::endl;
  }

  if (commitSuccessSpinDB)
  {
    std::cout << "Commit to SpinDB : SUCCESS" << std::endl;
  }
  else
  {
    std::cout << "Commit to SpinDB : FAILURE" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int XingShiftCal::Calibrate(const int final)
{
  CalculateCrossingShift(xingshift, scalercounts, success);
  if (!success)
  {
    std::cout << "CROSSING CALIBRATION FAILED." << std::endl;
    if (final)
    {
      std::cout << "DONE CALIBRATING." << std::endl;
    }
    return 0;
  }

  if (final)
  {
    // ostringstream comment;
    std::cout << "CROSSING CALIBRATION SUCCESS. XINGSHIFT: " << xingshift << std::endl;
    std::cout << "DONE CALIBRATING." << std::endl;
    // AddComment(comment.str());
  }

  return 0;
}

int XingShiftCal::CalculateCrossingShift(int &xing, uint64_t counts[NTRIG][NBUNCHES], bool &succ)
{
  succ = false;
  int shift_array[NTRIG] = {0};

  int trig_inactive_array[NTRIG] = {0};

  int last_active_index = 0;

  int _temp;
  for (int itrig = 0; itrig < NTRIG; itrig++)
  {
    long long _counts = 0;
    for (int ii = 0; ii < NBUNCHES; ii++)
    {
      _counts += counts[itrig][ii];
    }

    if (_counts < 10000)
    {
      trig_inactive_array[itrig] = 1;
    }
    else
    {
      last_active_index = itrig;
    }

    long long abort_sum_prev = _counts;

    _temp = 0;
    for (int ishift = 0; ishift < NBUNCHES; ishift++)
    {
      long long abort_sum = 0;
      for (int iunfillbunch = 0; iunfillbunch < NBUNCHES; iunfillbunch++)
      {
	if (blueFillPattern[iunfillbunch] && yellFillPattern[iunfillbunch])
	{
	  continue;
	}
	int shiftbunch = iunfillbunch - ishift;
	if (shiftbunch < 0)
	{
	  shiftbunch = 120 + shiftbunch;
	}
        abort_sum += counts[itrig][(shiftbunch) % NBUNCHES];
      }
      if (abort_sum < abort_sum_prev)
      {
        abort_sum_prev = abort_sum;
        _temp = ishift;
      }
    }

    shift_array[itrig] = _temp;
  }

  for (int itrig = 0; itrig < NTRIG; itrig++)
  {
    // if not matching for all trigger selections used, fails
    if (!trig_inactive_array[itrig])
    {
      if (shift_array[itrig] == shift_array[last_active_index])
      {
        xing = shift_array[itrig];
        succ = true;
      }
      else
      {
        xing = -999;
        succ = false;
        return 0;
      }
    }
  }

  // succ = true;
  return 0;
}

int XingShiftCal::WriteToCDB(const std::string &fname)
{
  std::cout << "XingShiftCal::WriteToCDB()" << std::endl;
  int xing_correction_offset = 0;
  if (success)
  {
    xing_correction_offset = xingshift;
    CDBTTree *cdbttree = new CDBTTree(fname);
    cdbttree->SetSingleIntValue("crossingshift", xing_correction_offset);
    cdbttree->CommitSingle();
    // cdbttree->Print();
    cdbttree->WriteCDBTTree();
    delete cdbttree;

    commitSuccessCDB = 1;
  }
  else
  {
    // if (verbosity) {
    std::cout << "no successful calibration, do not commit crossing shift to CDB" << std::endl;
    //}
    commitSuccessCDB = 0;
  }

  return 0;
}

// Commit to spinDB ported from PHENIX
int XingShiftCal::CommitToSpinDB()
{
  std::cout << "XingShiftCal::CommitPatternToSpinDB()" << std::endl;
  std::string status;  //-------------------------------------------->


  if (runnumber == 0)
  {
    std::cout << "Run doesn't exist" << std::endl;
    commitSuccessSpinDB = 0;
    return 0;
  }

  if (fillnumberBlue != fillnumberYellow)
  {
    std::cout << "fillnumber is wrong : fillnumberBlue = " << fillnumberBlue
              << "fillnumberYellow = " << fillnumberYellow << std::endl;
    commitSuccessSpinDB = 0;
    return 0;
  }

  

  int xing_correction_offset = -999;
  if (success)
  {
    xing_correction_offset = xingshift;
  }
  else
  {
    // if (verbosity) {
    std::cout << "no successful calibration, commit crossing shift -999 to spinDB" << std::endl;
    
    //}
    commitSuccessSpinDB = 0;
    //return 0;
  }

  // prepare values for db
  unsigned int qa_level = 0xffff;
  

  

  // if (verbosity) {
  std::cout << "polb = " << polBlue << " +- " << polBlueErr
            << std::endl
            << "poly = " << polYellow << " +- " << polYellowErr
            << std::endl;
  //}

  if (true /*0 && verbosity*/)
  {
    for (int ibeam = 0; ibeam < 2; ibeam++)
    {
      if (!ibeam)
      {
        std::cout << "spinpatternblue = {";
      }
      else
      {
        std::cout << "spinpatternyellow = {";
      }
      for (int icross = 0; icross < NBUNCHES; icross++)
      {
        if (!ibeam)
        {
          std::cout << blueSpinPattern[icross] << ",";
        }
        else
        {
          std::cout << yellSpinPattern[icross] << ",";
        }
      }
      std::cout << "\b}\n";
    }
  }

  // connect to spin db
  std::string dbname = "spinDB_write";
  std::string dbowner = "phnxrc";
  std::string dbpasswd = "";
  std::string dbtable = "spin";
  odbc::Connection *conSpin = nullptr;
  try
  {
    conSpin = odbc::DriverManager::getConnection(dbname.c_str(), dbowner.c_str(), dbpasswd.c_str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE
              << " Exception caught at XingShiftCal::CommitPatternToSpinDB when connecting to spin DB" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    commitSuccessSpinDB = 0;
    if (conSpin)
    {
      delete conSpin;
      conSpin = nullptr;
    }
    return 0;
  }
  // if (verbosity) cout << "opened spin DB connection" << endl;

  // check if this run already exists in spin_oncal
  bool runExists = false;
  std::ostringstream sqlSpinSelect;
  sqlSpinSelect << "SELECT runnumber, qa_level FROM " << dbtable
                << " WHERE runnumber = " << runnumber
                << " AND qa_level = " << qa_level
                << ";";
  // if (verbosity) cout<<sqlSpinSelect.str()<<endl;
  odbc::Statement *stmtSpinSelect = conSpin->createStatement();
  odbc::ResultSet *rsSpin = nullptr;
  try
  {
    rsSpin = stmtSpinSelect->executeQuery(sqlSpinSelect.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE
              << " Exception caught at XingShiftCal::CommitPatternToSpinDB when querying spin DB" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    commitSuccessSpinDB = 0;
    if (conSpin)
    {
      delete conSpin;
      conSpin = nullptr;
    }
    return 0;
  }
  if (rsSpin->next())
  {
    if (true /*verbosity*/)
    {
      std::cout << "run " << runnumber << " exists in " << dbtable
                << ", ready to UPDATE" << std::endl;
    }
    runExists = true;
  }
  else
  {
    if (true /*verbosity*/)
    {
      std::cout << "run " << runnumber << " NOT exists in " << dbtable
                << ", ready to INSERT" << std::endl;
    }
  }

  if (runExists && !overwriteSpinEntry)
  {
    std::cout << "BUT overwriteSpinEntry = " << overwriteSpinEntry << std::endl;
    std::cout << "XingShiftCal is NOT going to UPDATE the entry" << std::endl;
    commitSuccessSpinDB = 0;
    if (conSpin)
    {
      delete conSpin;
      conSpin = nullptr;
    }
    return 0;
  }

  // prepare insert sql
  std::ostringstream sql;
  if (runExists)
  {
    // SQLArrayConstF(float x, int n)
    sql << "UPDATE " << dbtable
        << " SET fillnumber = " << fillnumberBlue << ", "
        << "polarblue = " << SQLArrayConstF(polBlue, NBUNCHES) << ", "
        << "polarblueerror = " << SQLArrayConstF(polBlueErr, NBUNCHES) << ", "
        << "polaryellow = " << SQLArrayConstF(polYellow, NBUNCHES) << ", "
        << "polaryellowerror = " << SQLArrayConstF(polYellowErr, NBUNCHES) << ", "
	<< "crossingshift = " << xing_correction_offset << ", ";
    sql << "spinpatternblue = '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << blueSpinPattern[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";
    sql << ", spinpatternyellow = '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << yellSpinPattern[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";

    sql << " WHERE runnumber = " << runnumber
        << " AND qa_level = " << qa_level
        << ";";
  }
  else
  {
    sql << "INSERT INTO " << dbtable;
    sql << " (runnumber, fillnumber, polarblue, polarblueerror, polaryellow, polaryellowerror, crossingshift, spinpatternblue, spinpatternyellow, qa_level) VALUES (";
    sql << runnumber << ", "
        << fillnumberBlue << ", "
        << SQLArrayConstF(polBlue, NBUNCHES) << ", "
        << SQLArrayConstF(polBlueErr, NBUNCHES) << ", "
        << SQLArrayConstF(polYellow, NBUNCHES) << ", "
        << SQLArrayConstF(polYellowErr, NBUNCHES) << ", "
	<< xing_correction_offset << ", ";  
    sql << "'{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << blueSpinPattern[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";
    sql << ", '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << yellSpinPattern[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";
    sql << ", " << qa_level << ");";
  }
  if (true /*verbosity*/)
  {
    std::cout << sql.str() << std::endl;
  }

  // exec sql

  odbc::Statement *stmtSpin = conSpin->createStatement();
  try
  {
    stmtSpin->executeUpdate(sql.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE
              << " Exception caught at XingShiftCal::CommitPatternToSpinDB when insert into spin DB" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    commitSuccessSpinDB = 0;
    if (conSpin)
    {
      delete conSpin;
      conSpin = nullptr;
    }
    return 0;
  }

  // if (verbosity) cout<<"spin db done"<<endl;
  commitSuccessSpinDB = 1;

  if (conSpin)
  {
    delete conSpin;
    conSpin = nullptr;
  }
  return 0;
}

std::string XingShiftCal::SQLArrayConstF(float x, int n)
{
  std::ostringstream s;
  s << "'{";
  for (int i = 0; i < n; i++)
  {
    s << x;
    if (i < n - 1)
    {
      s << ",";
    }
  }
  s << "}'";
  return s.str();
}

//____________________________________________________________________________..
int XingShiftCal::Reset(PHCompositeNode * /*topNode*/)
{
  // std::cout << "XingShiftCal::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void XingShiftCal::Print(const std::string &what) const
{
  std::cout << "XingShiftCal::Print(const std::string &what) const Printing info for " << what << std::endl;
}
