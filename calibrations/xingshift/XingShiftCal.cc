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

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>
#include <odbc++/types.h>

#include <format>
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
  preset_pattern_blue["111x111_P1"] = "+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+*********";
  preset_pattern_yellow["111x111_P1"] = "++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++-*********";
  preset_pattern_blue["111x111_P2"] = "-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-*********";
  preset_pattern_yellow["111x111_P2"] = "++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++-*********";
  preset_pattern_blue["111x111_P3"] = "+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+*********";
  preset_pattern_yellow["111x111_P3"] = "--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--+*********";
  preset_pattern_blue["111x111_P4"] = "-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-*********";
  preset_pattern_yellow["111x111_P4"] = "--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--+*********";
  preset_pattern_blue["111x111_P5"] = "++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++-*********";
  preset_pattern_yellow["111x111_P5"] = "+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+*********";
  preset_pattern_blue["111x111_P6"] = "--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--+*********";
  preset_pattern_yellow["111x111_P6"] = "+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+--+-++-+-+-+--+-+-+-++-+*********";
  preset_pattern_blue["111x111_P7"] = "++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++-*********";
  preset_pattern_yellow["111x111_P7"] = "-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-*********";
  preset_pattern_blue["111x111_P8"] = "--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--++--+*********";
  preset_pattern_yellow["111x111_P8"] = "-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-++-+--+-+-+-++-+-+-+--+-*********";

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

    pBlueAsym = evt->getPacket(packet_BLUEASYM);
    pYellAsym = evt->getPacket(packet_YELLASYM);

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


     //============== Get pC spin patterns from buckets ==============//

    // Get bunch asymmetries for measured spin pattern
    // there are 360 buckets for 120 bunches
    if (pBlueAsym)
    {
      for (int i = 0; i < 360; i += 3)
      {
        float blueAsyms = pBlueAsym->iValue(i) / 10000.0;
        float blueAsymsErr = pBlueAsym->iValue(i + 360) / 10000.0;

        float bluebot = blueAsyms - blueAsymsErr;
        float bluetop = blueAsyms + blueAsymsErr;

        if (blueAsyms != 0 || bluebot != 0 || bluetop != 0)
        {
          if (bluebot > 0 && bluetop > 0)
          {
	    bluePcSpinPattern[i / 3] = 1;
          }
          else if (bluebot < 0 && bluetop < 0)
          {
            bluePcSpinPattern[i / 3] = -1;
          }
          else if (bluebot <= 0 && bluetop >= 0)
          {
            bluePcSpinPattern[i / 3] = 0;
          }
        }
        else
        {
          bluePcSpinPattern[i / 3] = 10;
        }
      }
      delete pBlueAsym;
    }

    if (pYellAsym)
    {
      for (int i = 0; i < 360; i += 3)
      {
        float yellAsyms = pYellAsym->iValue(i) / 10000.0;
        float yellAsymsErr = pYellAsym->iValue(i + 360) / 10000.0;

        float yellbot = yellAsyms - yellAsymsErr;
        float yelltop = yellAsyms + yellAsymsErr;

        if (yellAsyms != 0 || yellbot != 0 || yelltop != 0)
        {
          if (yellbot > 0 && yelltop > 0)
          {
            yellPcSpinPattern[i / 3] = 1;
          }
          else if (yellbot < 0 && yelltop < 0)
          {
            yellPcSpinPattern[i / 3] = -1;
          }
          else if (yellbot <= 0 && yelltop >= 0)
          {
            yellPcSpinPattern[i / 3] = 0;
          }
        }
        else
        {
          yellPcSpinPattern[i / 3] = 10;
        }
      }
      delete pYellAsym;
    }
    //=======================================================================//


    //============== Get fill number ==============//
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

  const std::string cdbfname = std::format("SPIN-{:08}_crossingshiftCDBTTree.root", runnumber);
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


  //=============== connect to daq db to get gl1p scalers ===============
  std::string daqdbname = "daq";
  std::string daqdbowner = "phnxrc";
  std::string daqdbpasswd = "";

  odbc::Connection *conDAQ = nullptr;
  try
  {
    conDAQ = odbc::DriverManager::getConnection(daqdbname.c_str(), daqdbowner.c_str(), daqdbpasswd.c_str());
  }
  catch (odbc::SQLException &eDAQ)
  {
    std::cout << PHWHERE
              << " Exception caught at XingShiftCal::CommitPatternToSpinDB when connecting to DAQ DB" << std::endl;
    std::cout << "Message: " << eDAQ.getMessage() << std::endl;
    //commitSuccessDAQDB = 0;
    if (conDAQ)
    {
      delete conDAQ;
      conDAQ = nullptr;
    }
    return 0;
  }

  std::ostringstream sqlGL1PSelect;
  sqlGL1PSelect << "SELECT index, bunch, scaled FROM gl1_pscalers"
                << " WHERE runnumber = " << runnumber
                << ";";
  odbc::Statement *stmtGL1PSelect = conDAQ->createStatement();
  odbc::ResultSet *rsGL1P = nullptr;
  try
  {
    rsGL1P = stmtGL1PSelect->executeQuery(sqlGL1PSelect.str());
  }
  catch (odbc::SQLException &eGL1P)
  {
    std::cout << PHWHERE
              << " Exception caught at XingShiftCal::CommitPatternToSpinDB when querying DAQ DB" << std::endl;
    std::cout << "Message: " << eGL1P.getMessage() << std::endl;
    //commitSuccessSpinDB = 0;
    if (conDAQ)
    {
      delete conDAQ;
      conDAQ = nullptr;
    }
    return 0;
  }


  while (rsGL1P->next()) {
    int index = rsGL1P->getInt("index");
    int bunch = rsGL1P->getInt("bunch");
    //MBD NS
    if (index == 0)
    {
      mbdns[bunch] = rsGL1P->getInt("scaled");
    }
    else if (index == 1)
    {
      mbdvtx[bunch] = rsGL1P->getInt("scaled");
    }
    else if (index == 5)
    {
      zdcns[bunch] = rsGL1P->getInt("scaled");
    }

  }
  // =======================================================
  

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

  if (true)
  {
    for (int i = 0; i < 3; i++)
    {
      if (i == 0)
      {
	std::cout << "mbdns = {";
      }
      else if (i == 1)
      {
	std::cout << "mbdvtx = {";
      }
      else if (i == 2)
      {
	std::cout << "zdcns = {";
      }
      for (int icross = 0; icross < NBUNCHES; icross++)
      {
        if (i == 0)
        {
          std::cout << mbdns[icross] << ",";
        }
        else if (i == 1)
        {
          std::cout << mbdvtx[icross] << ",";
        }
	else if (i == 2)
	{
	  std::cout << zdcns[icross] << ",";
	}	
      }
      std::cout << "\b}\n";
    }  
  }



  //================ connect to spin db ====================
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

    sql << ", mbdns = '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << mbdns[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";

    sql << ", mbdvtx = '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << mbdvtx[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";

    sql << ", zdcns = '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << zdcns[icross];
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
    sql << " (runnumber, fillnumber, polarblue, polarblueerror, polaryellow, polaryellowerror, crossingshift, spinpatternblue, spinpatternyellow, mbdns, mbdvtx, zdcns, qa_level) VALUES (";
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
    sql << ", '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << mbdns[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";
    sql << ", '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << mbdvtx[icross];
      if (icross < NBUNCHES - 1)
      {
        sql << ",";
      }
    }
    sql << "}'";
    sql << ", '{";
    for (int icross = 0; icross < NBUNCHES; icross++)
    {
      sql << zdcns[icross];
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

  if (conDAQ)
  {
    delete conDAQ;
    conDAQ = nullptr;
  }

  if (conSpin)
  {
    delete conSpin;
    conSpin = nullptr;
  }


  // ========== Do spin db qa here =========== //
  SpinDBQA();
  // ========================================= //


  return 0;
}


int XingShiftCal::SpinDBQA()
{

  // prepare values for db
  unsigned int qa_level = 0xffff;

  //================ connect to spin db ====================
  std::string dbname = "spinDB";
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

  // check if this run already exists in spin_oncal and get badrun value if it exists
  bool runExists = false;
  int prevbadrunval = 0;
  std::ostringstream sqlSpinSelect;
  sqlSpinSelect << "SELECT runnumber, qa_level, badrunqa FROM " << dbtable
                << " WHERE runnumber = " << runnumber
                << " AND qa_level = " << qa_level
                << ";";

  std::cout << "SELECT runnumber, qa_level, badrunqa FROM " << dbtable
                << " WHERE runnumber = " << runnumber
                << " AND qa_level = " << qa_level
	    << ";" << std::endl;
  
  odbc::Statement *stmtSpinSelect = conSpin->createStatement();
  odbc::ResultSet *rsSpin = nullptr;
  try
  {
    rsSpin = stmtSpinSelect->executeQuery(sqlSpinSelect.str());
  }
  catch (odbc::SQLException &e)
  {
    std::cout << PHWHERE
              << " Exception caught at XingShiftCal::SpinDBQA when querying spin DB" << std::endl;
    std::cout << "Message: " << e.getMessage() << std::endl;
    if (conSpin)
    {
      delete conSpin;
      conSpin = nullptr;
    }
    return 0;
  }
  if (rsSpin->next())
  {
    prevbadrunval = rsSpin->getInt("badrunqa");
    if (rsSpin->wasNull()) 
    {
      std::cout << "SPINDBQA: badrunqa is NULL. Setting to 0 for qa check. "<< std::endl;
      prevbadrunval = 0;
    }
    if (true)
    {
      std::cout << "run " << runnumber << " exists in " << dbtable
                << ", ready to do QA" << std::endl;
    }
    runExists = true;
  }
  else
  {
    if (true)
    {
      std::cout << "run " << runnumber << " NOT exists in " << dbtable
                << ", no QA" << std::endl;
    }
  }

  int badrunQA = 0;
  
  if (prevbadrunval > 0)
  {
    std::cout << "SPINDBQA: badrunqa is already > 0. No additional qa is performed." << std:: endl;
    if (conSpin)
    {
      delete conSpin;
      conSpin = nullptr;
    }
    return 0;  
  }
  else
  {
    // =========== Do bad run QA here =============
    // if (conditions pass && previously existing badRunQA != 1): badrunQA = 0
    // if (conditions fail (bad run)): badrunQA = 1

    //if spin pattern does not match known MCR pattern
    std::string scdev_blue = "";
    std::string scdev_yell = "";
    for (int crossing = 0; crossing < 120; crossing++)
    {
      int ibluespin = blueSpinPattern[crossing];
      int iyellspin = yellSpinPattern[crossing];
      if (ibluespin == 1)
      {
	scdev_blue.push_back('+');
      }
      else if (ibluespin == -1)
      {
	scdev_blue.push_back('-');
      }
      else
      {
	scdev_blue.push_back('*');
      }

      if (iyellspin == 1)
      {
	scdev_yell.push_back('+');
      }
      else if (iyellspin == -1)
      {
	scdev_yell.push_back('-');
      }
      else
      {
	scdev_yell.push_back('*');
      }
    }
    //std::string scdev_blue = TH1_to_string(hspinpatternBlue);
    //std::string scdev_yell = TH1_to_string(hspinpatternYellow);
    std::string pattern_name = "UNKNOWN";
    
    for (std::map<std::string, std::string>::const_iterator ii = preset_pattern_blue.begin(); ii != preset_pattern_blue.end(); ++ii)
    {
      std::string key = (*ii).first;
      if (preset_pattern_blue[key] == scdev_blue && preset_pattern_yellow[key] == scdev_yell)
      {
	pattern_name = key;
      }
    }

    if (pattern_name == "UNKNOWN")
    {
      badrunQA = 1;
      std::cout << "SPINDBQA: Pattern is unidentified from known CDEV pattern. Setting bad run." << std::endl;
    }


    //if pc spin pattern does not match intended spin pattern within < 10 bunches
    int mismatches = 0;
//    std::string spin_pattern_blue = "";
//    std::string spin_pattern_yell = "";
    for (int crossing = 0; crossing < 120; crossing++)
    {
      int spin_cdev_blue = blueSpinPattern[crossing];
      int spin_cdev_yell = yellSpinPattern[crossing];
      int spin_pC_blue = bluePcSpinPattern[crossing];
      int spin_pC_yell = yellPcSpinPattern[crossing];

      if(spin_pC_blue==-1 || spin_pC_blue==1)
      {
	if (spin_cdev_blue != spin_pC_blue && !(spin_cdev_blue == 0 && spin_pC_blue == 10))
	{
	  mismatches += 1;
	}
      }

      if(spin_pC_yell==-1 || spin_pC_yell==1)
      {
	if (spin_cdev_yell != spin_pC_yell && !(spin_cdev_blue == 0 && spin_pC_blue == 10))
	{
	  mismatches += 1;
	}
      } 
    }

    if (mismatches > 10)
    {
      badrunQA = 1;
      std::cout << "SPINDBQA: CDEV pattern has > 10 mismatched bunches from pC polarimeter. Setting bad run." << std::endl;
    }

    //if crossing shift != 0
    int xing_correction_offset = -999;
    if (success)
    {
      xing_correction_offset = xingshift;
    }

    if (xing_correction_offset != 0)
    {
      badrunQA = 1;
      std::cout << "SPINDBQA: Crossing shift does not equal 0. Setting bad run." << std::endl;
    }

    //if polarization <= 0 || > 1.00 //makes sure polarization from CNI aren't garbage values
    if (polBlue <= 0.0 || polBlue > 100.0)
    {
      badrunQA = 1;
      std::cout << "SPINDBQA: Blue beam polarization is unknown. Setting bad run." << std::endl;
    }
    if (polYellow <= 0.0 || polYellow > 100.0)
    {
      badrunQA = 1;
      std::cout << "SPINDBQA: Yellow beam polarization is unknown. Setting bad run." << std::endl;
    }
      
    // ============================================
    
  }
  

  //================ connect to spin db write ====================
  dbname = "spinDB_write";
  dbowner = "phnxrc";
  dbpasswd = "";
  dbtable = "spin";
  conSpin = nullptr;
  if (runExists)
  {

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


    std::cout << "UPDATE " << dbtable
        << " SET badrunqa = " << badrunQA
	<< " WHERE runnumber = " << runnumber
        << " AND qa_level = " << qa_level
	      << ";" << std::endl;

    // prepare insert sql
    std::ostringstream sql;
    sql << "UPDATE " << dbtable
        << " SET badrunqa = " << badrunQA
	<< " WHERE runnumber = " << runnumber
        << " AND qa_level = " << qa_level
        << ";";

    // exec sql
    odbc::Statement *stmtSpin = conSpin->createStatement();
    try
    {
      stmtSpin->executeUpdate(sql.str());
    }
    catch (odbc::SQLException &e)
    {
      std::cout << PHWHERE
		<< " Exception caught at XingShiftCal::SpinDBQA when insert badrunqa into spin DB" << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
      commitSuccessSpinDB = 0;
      if (conSpin)
      {
	delete conSpin;
	conSpin = nullptr;
      }
      return 0;
    }

  }

  
  

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

