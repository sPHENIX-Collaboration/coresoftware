#include "CDBInterface.h"

#include <sphenixnpc/SphenixClient.h>

#include <ffaobjects/CdbUrlSave.h>
#include <ffaobjects/CdbUrlSavev1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <cstdint>  // for uint64_t
#include <filesystem>
#include <fstream>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <sstream>
#include <utility>  // for pair
#include <vector>   // for vector

CDBInterface *CDBInterface::__instance{nullptr};

CDBInterface *CDBInterface::instance()
{
  if (__instance)
  {
    return __instance;
  }
  __instance = new CDBInterface();
  return __instance;
}

//____________________________________________________________________________..
CDBInterface::CDBInterface(const std::string &name)
  : SubsysReco(name)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->addNewSubsystem(this);
}

//____________________________________________________________________________..

CDBInterface::~CDBInterface()
{
  delete cdbclient;
}

//____________________________________________________________________________..
int CDBInterface::End(PHCompositeNode *topNode)
{
  int iret = UpdateRunNode(topNode);
  PHNodeIterator iter(topNode);
  return iret;
}

//____________________________________________________________________________..
int CDBInterface::UpdateRunNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  CdbUrlSave *cdburls = findNode::getClass<CdbUrlSave>(runNode, "CdbUrl");
  if (!cdburls)
  {
    cdburls = new CdbUrlSavev1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cdburls, "CdbUrl", "PHObject");
    runNode->addNode(newNode);
  }
  else
  {
    std::set<std::tuple<std::string, std::string, uint64_t>> tmp_set;
    for (const auto &cdburl : *cdburls)
    {
      tmp_set.insert(cdburl);
    }
    // remove duplicates in our set
    // not possible using for range loops, iterator gets invalidated
    for (auto itr = m_UrlVector.cbegin(); itr != m_UrlVector.cend();)
    {
      if (tmp_set.contains(*itr))
      {
        if (Verbosity() > 2)
        {
          std::cout << PHWHERE << " cleaning duplicately saved: domain " << std::get<0>(*itr)
                    << ", url: " << std::get<1>(*itr)
                    << ", timestamp: " << std::get<2>(*itr) << std::endl;
        }
        itr = m_UrlVector.erase(itr);
      }
      else
      {
        ++itr;
      }
    }
  }
  for (const auto &tuple : m_UrlVector)
  {
    cdburls->AddUrl(tuple);
  }
  if (Verbosity() > 1)
  {
    cdburls->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CDBInterface::Print(const std::string & /* what */) const
{
  for (const auto &iter : m_UrlVector)
  {
    std::cout << "domain: " << std::get<0>(iter)
              << ", url: " << std::get<1>(iter)
              << ", timestamp: " << std::get<2>(iter) << std::endl;
  }
}

std::string CDBInterface::getUrl(const std::string &domain, const std::string &filename)
{
  if (disable)
  {
    return "";
  }
  std::string domain_noconst = domain;
  if (m_Read_From_File_Flag)
  {
    if (m_Payload_Url_Cache.contains(domain_noconst))
    {
      return m_Payload_Url_Cache[domain_noconst];
    }
    std::cout << "calibration " << domain << " not found in local cache" << std::endl;
    return "";
  }
  recoConsts *rc = recoConsts::instance();
  if (!rc->FlagExist("CDB_GLOBALTAG"))
  {
    std::cout << PHWHERE << "CDB_GLOBALTAG flag needs to be set via" << std::endl;
    std::cout << "rc->set_StringFlag(\"CDB_GLOBALTAG\",<global tag>)" << std::endl;
    gSystem->Exit(1);
  }
  if (!rc->FlagExist("TIMESTAMP"))
  {
    std::cout << PHWHERE << "TIMESTAMP flag needs to be set via" << std::endl;
    std::cout << "rc->set_uint64Flag(\"TIMESTAMP\",<64 bit timestamp>)" << std::endl;
    gSystem->Exit(1);
  }
  if (cdbclient == nullptr)
  {
    cdbclient = new SphenixClient(rc->get_StringFlag("CDB_GLOBALTAG"));
  }
  uint64_t timestamp = rc->get_uint64Flag("TIMESTAMP");
  if (Verbosity() > 0)
  {
    std::cout << "Global Tag: " << rc->get_StringFlag("CDB_GLOBALTAG")
              << ", domain: " << domain_noconst
              << ", timestamp: " << timestamp;
  }
  std::string return_url = cdbclient->getCalibration(domain_noconst, timestamp);
  if (return_url.empty())
  {
    if (!disable_default)
    {
      std::string domain_copy = domain_noconst;
      domain_noconst = domain_noconst + "_default";
      return_url = cdbclient->getCalibration(domain_noconst, timestamp);
      if (return_url.empty())
      {
        if (Verbosity() > 0)
        {
          std::cout << "... reply: no file found for "
                    << domain_copy << " or " << domain_noconst << std::endl;
        }
        return_url = filename;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        std::cout << "... reply: no file found for "
                  << domain_noconst << std::endl;
      }
    }
  }
  if (Verbosity() > 0)
  {
    if (!return_url.empty())
    {
      std::cout << "... reply: " << return_url << std::endl;
    }
  }
  if (!return_url.empty())
  {
    auto pret = m_UrlVector.insert(make_tuple(domain_noconst, return_url, timestamp));
    if (!pret.second && Verbosity() > 1)
    {
      std::cout << PHWHERE << "not adding again " << domain_noconst << ", url: " << return_url
                << ", time stamp: " << timestamp << std::endl;
    }
  }
  return return_url;
}

void CDBInterface::DumpCalibrations(const std::string &filename)
{
  recoConsts *rc = recoConsts::instance();
  if (!rc->FlagExist("CDB_GLOBALTAG"))
  {
    std::cout << PHWHERE << "CDB_GLOBALTAG flag needs to be set via" << std::endl;
    std::cout << "rc->set_StringFlag(\"CDB_GLOBALTAG\",<global tag>)" << std::endl;
    gSystem->Exit(1);
  }
  if (!rc->FlagExist("TIMESTAMP"))
  {
    std::cout << PHWHERE << "TIMESTAMP flag needs to be set via" << std::endl;
    std::cout << "rc->set_uint64Flag(\"TIMESTAMP\",<64 bit timestamp>)" << std::endl;
    gSystem->Exit(1);
  }
  if (cdbclient == nullptr)
  {
    cdbclient = new SphenixClient(rc->get_StringFlag("CDB_GLOBALTAG"));
  }
  uint64_t timestamp = rc->get_uint64Flag("TIMESTAMP");
  cdbclient->DumpCalibrations(timestamp, filename);
  return;
}

void CDBInterface::ReadCalibrationsFromFile(const std::string &filename)
{
  std::filesystem::path filePath = filename;
  if (!std::filesystem::exists(filePath))
  {
    std::cout << PHWHERE << " cannot locate " << filename << std::endl;
    gSystem->Exit(1);
  }
  if (!std::filesystem::is_regular_file(filePath))
  {
    std::cout << PHWHERE << " not a regular file " << filename << std::endl;
    gSystem->Exit(1);
  }
  std::ifstream calibsfile(filename);
  if (calibsfile.is_open())
  {
    std::string line;
    while (std::getline(calibsfile, line))
    {
      // Skip empty lines
      if (line.empty())
      {
        continue;
      }

      // Skip comments
      if (line[0] == '#')
      {
        continue;
      }
      std::istringstream iss(line);
      std::string key;
      std::string payload_url;
      if (iss >> key >> payload_url)
      {
        m_Payload_Url_Cache.insert(std::make_pair(key, payload_url));
      }
    }
    m_Read_From_File_Flag = true;
  }
  else
  {
    std::cout << "could not open " << filename << std::endl;
  }
  return;
}
