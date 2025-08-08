#include "TriggerRunInfoReco.h"
#include "TriggerRunInfov1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <odbc++/connection.h>
#include <odbc++/drivermanager.h>
#include <odbc++/resultset.h>
#include <odbc++/statement.h>
#include <odbc++/types.h>

#include <array>
#include <cstdint>
#include <iostream>

TriggerRunInfoReco::TriggerRunInfoReco(const std::string &name)
  : SubsysReco(name)
{
}

int TriggerRunInfoReco::Init(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << "RUN node not found!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Create the TriggerRunInfo object and add it to the RUN node
  TriggerRunInfo *triggerRunInfo = new TriggerRunInfov1();
  PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(triggerRunInfo, "TriggerRunInfo", "PHObject");
  runNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TriggerRunInfoReco::InitRun(PHCompositeNode *topNode)
{
  recoConsts *rc = recoConsts::instance();
  int runnumber = rc->get_IntFlag("RUNNUMBER");

  // Retrieve the TriggerRunInfo object
  TriggerRunInfo *triggerRunInfo = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerRunInfo)
  {
    std::cout << "TriggerRunInfo object not found!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (m_useEmulator)
  {
    SetTriggerEmulator(triggerRunInfo);
    if (Verbosity())
    {
      triggerRunInfo->identify();
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }
  // Fetch trigger prescales and fill the TriggerRunInfo object
  if (fetchTriggerPrescales(runnumber, triggerRunInfo) != 0)
  {
    std::cout << "Failed to fetch trigger prescales for run number " << runnumber << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Fetch trigger scalers and fill the TriggerRunInfo object
  if (fetchTriggerScalers(runnumber, triggerRunInfo) != 0)
  {
    std::cout << "Failed to fetch trigger prescales for run number " << runnumber << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (Verbosity())
  {
    triggerRunInfo->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
void TriggerRunInfoReco::SetTriggerEmulator(TriggerRunInfo *triggerRunInfo)
{
  std::string names[64];
  for (int i = 0; i < 64; i++)
  {
    names[i] = "unknown" + std::to_string(i);
  }
  for (int i = 0; i < 64; i++)
  {
    triggerRunInfo->setTrigger(i, "unknown", i, 0);
  }

  triggerRunInfo->setTrigger(28, "Photon 2 GeV", 28, 1);
  triggerRunInfo->setTrigger(29, "Photon 3 GeV", 29, 1);
  triggerRunInfo->setTrigger(30, "Photon 4 GeV", 30, 1);
  triggerRunInfo->setTrigger(31, "Photon 5 GeV", 31, 1);
  triggerRunInfo->setTrigger(20, "Jet 6 GeV", 20, 1);
  triggerRunInfo->setTrigger(21, "Jet 8 GeV", 21, 1);
  triggerRunInfo->setTrigger(22, "Jet 10 GeV", 22, 1);
  triggerRunInfo->setTrigger(23, "Jet 12 GeV", 23, 1);

  return;
}

int TriggerRunInfoReco::fetchTriggerPrescales(int runnumber, TriggerRunInfo *triggerRunInfo)
{
  odbc::Connection *dbConnection = nullptr;
  try
  {
    dbConnection = odbc::DriverManager::getConnection("daq", "", "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Database connection failed: " << e.getMessage() << std::endl;
    return 1;
  }

  std::string sql = "SELECT * FROM gl1_scaledown WHERE runnumber = " + std::to_string(runnumber) + ";";
  odbc::Statement *stmt = dbConnection->createStatement();
  odbc::ResultSet *resultSet = stmt->executeQuery(sql);
  int prescales[64];
  std::string names[64];
  for (int i = 0; i < 64; i++)
  {
    names[i] = "unknown" + std::to_string(i);
  }

  if (resultSet && resultSet->next())
  {
    // Iterate over the columns and fill the TriggerRunInfo object
    for (int bit = 0; bit < 64; ++bit)
    {
      std::string columnName = std::string("scaledown") + (bit < 10 ? "0" : "") + std::to_string(bit);
      prescales[bit] = resultSet->getInt(columnName);
    }
  }
  else
  {
    std::cout << "No data found for run number " << runnumber << std::endl;
    delete resultSet;
    delete stmt;
    delete dbConnection;
    return 1;
  }

  sql = "SELECT * FROM gl1_triggernames WHERE runnumber_last > " + std::to_string(runnumber) + " AND runnumber < " + std::to_string(runnumber) + ";";
  stmt = dbConnection->createStatement();
  resultSet = stmt->executeQuery(sql);

  while (resultSet && resultSet->next())
  {
    int bit = resultSet->getInt("index");
    names[bit] = resultSet->getString("triggername");
  }

  for (int bit = 0; bit < 64; bit++)
  {
    triggerRunInfo->setTrigger(bit, names[bit], bit, prescales[bit]);
  }

  delete resultSet;
  delete stmt;
  delete dbConnection;
  return 0;
}

int TriggerRunInfoReco::fetchTriggerScalers(int runnumber, TriggerRunInfo *triggerRunInfo)
{
  odbc::Connection *dbConnection = nullptr;
  try
  {
    dbConnection = odbc::DriverManager::getConnection("daq", "", "");
  }
  catch (odbc::SQLException &e)
  {
    std::cout << "Database connection failed: " << e.getMessage() << std::endl;
    return 1;
  }

  std::string sql = "SELECT * FROM gl1_scalers WHERE runnumber = " + std::to_string(runnumber) + ";";
  odbc::Statement *stmt = dbConnection->createStatement();
  odbc::ResultSet *resultSet = stmt->executeQuery(sql);
  std::array<std::array<uint64_t, 3>, 64> scalers{};  // initialize to zero
  if (!resultSet)
  {
    std::cout << "No data found for run number " << runnumber << std::endl;
    delete resultSet;
    delete stmt;
    delete dbConnection;
    return 1;
  }

  while (resultSet->next())
  {
    int index = resultSet->getInt("index");
    // Iterate over the columns and fill the TriggerRunInfo object
    scalers[index][0] = resultSet->getLong("scaled");
    scalers[index][1] = resultSet->getLong("live");
    scalers[index][2] = resultSet->getLong("raw");
  }

  delete resultSet;
  delete stmt;
  delete dbConnection;

  for (int i = 0; i < 64; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      triggerRunInfo->setTriggerScalers(i, j, scalers[i][j]);
    }
    double scaled = static_cast<double>(scalers[i][0]);
    double live = static_cast<double>(scalers[i][1]);
    double prescale = -1;
    if (scaled >= 1)
    {
      prescale = live / scaled;
    }
    triggerRunInfo->setTriggerPrescale(i, prescale);
  }

  return 0;
}
