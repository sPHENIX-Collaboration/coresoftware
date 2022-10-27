#include "XploadInterface.h"

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
#include <phool/recoConsts.h>

#include <xpload/xpload.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair
#include <vector>    // for vector

XploadInterface *XploadInterface::__instance = nullptr;

XploadInterface *XploadInterface::instance()
{
  if (__instance)
  {
    return __instance;
  }
  __instance = new XploadInterface();
  return __instance;
}

//____________________________________________________________________________..
XploadInterface::XploadInterface(const std::string &name)
  : SubsysReco(name)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->addNewSubsystem(this);
}

//____________________________________________________________________________..
int XploadInterface::End(PHCompositeNode *topNode)
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
    for (const auto & cdburl : *cdburls)
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
          std::cout << "removing already saved: domain " << std::get<0>(*itr)
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
  if (Verbosity() > 0)
  {
    cdburls->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void XploadInterface::Print(const std::string & /* what */) const
{
  for (auto &iter : m_UrlVector)
  {
    std::cout << "domain: " << std::get<0>(iter)
              << ", url: " << std::get<1>(iter)
              << ", timestamp: " << std::get<2>(iter) << std::endl;
  }
}

std::string XploadInterface::getUrl(const std::string &domain, const std::string &filename)
{
  std::string return_url = filename;
  recoConsts *rc = recoConsts::instance();
  uint64_t timestamp = rc->get_uint64Flag("TIMESTAMP", 12345678912345);
  xpload::Result result = xpload::fetch(rc->get_StringFlag("XPLOAD_TAG", "TEST"), domain, timestamp, xpload::Configurator(rc->get_StringFlag("XPLOAD_CONFIG", "sPHENIX_cdb")));
  if (!result.payload.empty())
  {
    return_url = result.payload;
  }
  auto pret = m_UrlVector.insert(make_tuple(domain, return_url, timestamp));
  if (!pret.second && Verbosity() > 1)
  {
    std::cout << "duplicate entry " << domain << ", url: " << return_url
              << ", time stamp: " << timestamp << std::endl;
  }
  return return_url;
}
