#include "CaloTemp.h"

#include <boost/format.hpp>

// Tower includes
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TProfile2D.h>

#include <Event/Event.h>
#include <Event/packet.h>
#include <cassert>
#include <sstream>
#include <string>
#include <chrono>
#include <iostream>
#include <thread>

#include <ffaobjects/RunHeader.h>

#include <odbc++/connection.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <odbc++/drivermanager.h>
#pragma GCC diagnostic pop
#include <odbc++/resultset.h>
#include <odbc++/statement.h>  // for Statement
#include <odbc++/types.h>

using namespace std;

CaloTemp::CaloTemp(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(filename)
{
}

CaloTemp::~CaloTemp()
{
  delete hm;
}



int CaloTemp::Init(PHCompositeNode* /*topNode*/)
{
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here

  if (Verbosity() > 1) {
    std::cout << "Enter CaloTemp Init function" << std::endl;
  }

  outfile = new TFile(outfilename.c_str(), "UPDATE");

  if (detector == "CEMC") {
    h_calo_temp = new TProfile2D("h_cemc_sipm_temp", ";eta;phi", 96, 0, 96, 256, 0, 256);
    h_calo_temp2 = new TProfile2D("h_cemc_pa_temp", ";eta;phi", 96, 0, 96, 256, 0, 256);
  } else if (detector == "HCALIN") {
    h_calo_temp = new TProfile2D("h_ihcal_temp", ";eta;phi", 24, 0, 24, 64, 0, 64);
  } else if (detector == "HCALOUT") {
    h_calo_temp = new TProfile2D("h_ohcal_temp", ";eta;phi", 24, 0, 24, 64, 0, 64);
  } else {
    std::cout << "Unknown detector type " << detector << ". Options are CEMC, HCALIN and HCALOUT." << std::endl;
    return 1;
  }

  return 0;
}

int CaloTemp::InitRun(PHCompositeNode* topNode) 
{
  if (Verbosity() > 1) { 
    std::cout << "Get run header" << std::endl;
  }
  
  runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!runheader) {
    std::cout << "can't find runheader" << std::endl;
    return 1;
  }
  runnumber = runheader->get_RunNumber();

  if (Verbosity()) {
    std::cout << "run number " << runnumber << std::endl;
  }
  
  int result = getRunTime();
  if (result == 1) {
    std::cout << "Unable to connect to database. Exiting now." << std::endl;
    return 1;
  }
  
  result = getTempHist();
  if (result == 1) {
    std::cout << "Unable to fill calorimeter temperature from database. Exiting now." << std::endl;
    return 1;
  }

  if (Verbosity()) {
    std::cout << "runtime " << runtime << std::endl;
  }

  return 0;
}

int CaloTemp::getRunTime() {
  odbc::Connection *m_OdbcConnection = nullptr;
  int icount = 0;
  do {
    try {
      m_OdbcConnection = odbc::DriverManager::getConnection("daq","","");
    } catch (odbc::SQLException &e) {
      std::cout << " Exception caught during DriverManager::getConnection" << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
    }
    icount++;
    if (!m_OdbcConnection) { std::this_thread::sleep_for(std::chrono::seconds(30)); }  // sleep 30 seconds before retry
  } while (!m_OdbcConnection && icount < 3);
  if (icount == 3) { return 1; }

  std::string sql = "SELECT brtimestamp FROM run WHERE runnumber = " + std::to_string(runnumber) + ";";
  
  if (Verbosity() > 1) {
    std::cout << sql << std::endl;
  }
  odbc::Statement* stmt = m_OdbcConnection->createStatement();
  odbc::ResultSet* resultSet = stmt->executeQuery(sql);

  if (resultSet && resultSet->next()) {
    odbc::Timestamp brtimestamp = resultSet->getTimestamp("brtimestamp");
    runtime = brtimestamp.toString(); // brtimestamp is in 'America/New_York' time zone
  }

  if (Verbosity() > 1) { 
    std::cout << "runtime " << runtime << std::endl;
  }

  delete resultSet; 
  delete stmt; 
  delete m_OdbcConnection;
  return 0;

}


int CaloTemp::getTempHist() {
  
  odbc::Connection *m_OdbcConnection = nullptr;
  int icount = 0;
  do {
    try {
      m_OdbcConnection = odbc::DriverManager::getConnection("daq","","");
    } catch (odbc::SQLException &e) {
      std::cout << " Exception caught during DriverManager::getConnection" << std::endl;
      std::cout << "Message: " << e.getMessage() << std::endl;
    }
    icount++;
    if (!m_OdbcConnection) { std::this_thread::sleep_for(std::chrono::seconds(30)); }  // sleep 30 seconds before retry
  } while (!m_OdbcConnection && icount < 3);
  if (icount == 3) { return 1; }

  int det = 0;
  std::string tablename;
  std::string tempstring = "temp";
  if (detector == "CEMC") {
    tablename = "emcal_heartbeat";
    tempstring = "temp_sipm as temp, temp_pa";
  } else if (detector == "HCALIN") {
    tablename = "hcal_heartbeat";
    det = 1; 
  } else if (detector == "HCALOUT") {
    tablename = "hcal_heartbeat";
    det = 0;
  } else {
    std::cout << "Unknown detector type " << detector << ". Options are CEMC, HCALIN and HCALOUT." << std::endl;
    return 1;
  }

  // first get the closest time for the database 
  std::string sql = "SELECT time FROM " + tablename + " ORDER BY ABS(EXTRACT(epoch FROM time) - EXTRACT(epoch FROM '" + runtime + "'::timestamp)) LIMIT 1;";
  if (Verbosity() > 1) {
    std::cout << sql << std::endl;
  }

  odbc::Statement* timeStmt = m_OdbcConnection->createStatement();
  odbc::ResultSet* timeResultSet = timeStmt->executeQuery(sql);

  std::string closest_time;
  if (timeResultSet && timeResultSet->next()) {
    odbc::Timestamp timestamp = timeResultSet->getTimestamp("time");
    closest_time = timestamp.toString(); // Convert timestamp to string
  }

  if (Verbosity() > 1) {
    std::cout << "closest time " << closest_time << std::endl;
  }

  delete timeResultSet;
  delete timeStmt;

  odbc::Statement* tempStmt = m_OdbcConnection->createStatement();
  sql = "SELECT towerid, " + tempstring + " FROM " + tablename + " WHERE time = '" + closest_time + "'"; 
  if (detector == "HCALIN" || detector == "HCALOUT") {
    sql += " and detector = " + std::to_string(det);
  }
  sql += ";";
  
  if (Verbosity() > 1) {
    std::cout << sql << std::endl;
  }

  odbc::ResultSet* tempResultSet = tempStmt->executeQuery(sql);
  string calo_title = boost::str(boost::format("%s Temperature RunNumber %d RunTime %s DBTime %s") %detector %runnumber %runtime %closest_time);
  h_calo_temp->SetTitle(calo_title.c_str());
  if (detector == "CEMC") { 
    calo_title = boost::str(boost::format("%s SIPM Temperature RunNumber %d RunTime %s DBTime %s") %detector %runnumber %runtime %closest_time);
    string calo_title2 = calo_title = boost::str(boost::format("%s PreAmp Temperature RunNumber %d RunTime %s DBTime %s") %detector %runnumber %runtime %closest_time);
    h_calo_temp->SetTitle(calo_title.c_str());
    h_calo_temp2->SetTitle(calo_title2.c_str()); 
  }
  if (tempResultSet) {
    while (tempResultSet->next()) {
      int towerid = tempResultSet->getInt("towerid"); 
      float temp = tempResultSet->getFloat("temp"); 
      int calo_key;
      if (detector == "CEMC") { calo_key = TowerInfoDefs::encode_emcal(towerid); } 
      else { calo_key = TowerInfoDefs::encode_hcal(towerid); }
      int etabin = TowerInfoDefs::getCaloTowerEtaBin(calo_key);
      int phibin = TowerInfoDefs::getCaloTowerPhiBin(calo_key);
      h_calo_temp->Fill(etabin, phibin, temp);
      if (detector == "CEMC") { 
        float temp_pa = tempResultSet->getFloat("temp_pa");
        h_calo_temp2->Fill(etabin, phibin, temp_pa); 
      }
    }
} else {
    std::cerr << "Error: tempResultSet is NULL." << std::endl;
}

  delete tempResultSet;
  delete tempStmt;
  delete m_OdbcConnection;
  return 0;
}


int CaloTemp::End(PHCompositeNode* /*topNode*/)
{
  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}
