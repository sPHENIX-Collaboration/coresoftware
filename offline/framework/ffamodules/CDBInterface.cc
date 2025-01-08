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

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair
#include <vector>    // for vector

CDBInterface *CDBInterface::__instance = nullptr;

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
      if (tmp_set.find(*itr) != tmp_set.end())
      {
        if (Verbosity())
        {
          std::cout << PHWHERE << " removing already saved: domain " << std::get<0>(*itr)
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
  for (auto &tuple : m_UrlVector)
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
  for (auto &iter : m_UrlVector)
  {
    std::cout << "domain: " << std::get<0>(iter)
              << ", url: " << std::get<1>(iter)
              << ", timestamp: " << std::get<2>(iter) << std::endl;
  }
}

std::string CDBInterface::getUrl(const std::string &domain, const std::string &filename)
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
  if (Verbosity() > 0)
  {
    std::cout << "Global Tag: " << rc->get_StringFlag("CDB_GLOBALTAG")
              << ", domain: " << domain
              << ", timestamp: " << timestamp;
  }
  std::string return_url = cdbclient->getCalibration(domain, timestamp);
  if (Verbosity() > 0)
  {
    if (return_url.empty())
    {
      std::cout << "... reply: no file found" << std::endl;
    }
    else
    {
      std::cout << "... reply: " << return_url << std::endl;
    }
  }
  if (return_url.empty())
  {
    return_url = filename;
  }
  auto pret = m_UrlVector.insert(make_tuple(domain, return_url, timestamp));
  if (!pret.second && Verbosity() > 1)
  {
    std::cout << PHWHERE << "not adding again " << domain << ", url: " << return_url
              << ", time stamp: " << timestamp << std::endl;
  }
  return return_url;
}
