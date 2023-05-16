#include "CDBUtils.h"

#include "sphenixnpc.h"

#include <nlohmann/json.hpp>

#include <iostream>
#include <stdexcept>

CDBUtils::CDBUtils()
  : cdbclient(new sphenixnpc())
{}

CDBUtils::CDBUtils(const std::string &globaltag)
  : cdbclient(new sphenixnpc(globaltag))
{
}


void CDBUtils::createGlobalTag(const std::string &tagname)
{
  cdbclient->createGlobalTag(tagname);
}

int CDBUtils::deleteGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->deleteGlobalTag(tagname);
  int iret = resp["code"];
  if (iret != 0)
  {
    nlohmann::json msgcont = resp["msg"];
    std::cout << "message: " << msgcont << std::endl;
  }
  return iret;
}

void CDBUtils::lockGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->lockGlobalTag(tagname);
  int iret = resp["code"];
  if (iret != 0)
  {
    nlohmann::json msgcont = resp["msg"];
    std::cout << "message: " << msgcont << std::endl;
  }
  return;
}

void CDBUtils::unlockGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->unlockGlobalTag(tagname);
  int iret = resp["code"];
  if (iret != 0)
  {
    nlohmann::json msgcont = resp["msg"];
    std::cout << "message: " << msgcont << std::endl;
  }
  return;
}

void CDBUtils::clearCache()
{
  cdbclient->clearCache();
}

std::string CDBUtils::getUrl(const std::string &type, uint64_t iov)
{
  return cdbclient->getUrl(type, iov);
}

int CDBUtils::createDomain(const std::string &domain)
{
  cdbclient->createDomain(domain);
  int iret = -1;
  nlohmann::json resp;
  if (m_DomainCache.empty())
  {
    resp = cdbclient->getPayloadTypes();
    nlohmann::json msgcont = resp["msg"];
    for (auto &it : msgcont.items())
    {
      std::string existent_domain = it.value().at("name");
      m_DomainCache.insert(existent_domain);
    }
  }
  if (m_DomainCache.find(domain) == m_DomainCache.end())
  {
    resp = cdbclient->createPayloadType(domain);
    iret = resp["code"];
    if (iret == 0)
    {
      m_DomainCache.insert(domain);
    }
  }
  else
  {
    iret = 0;
  }
  return iret;
}

void CDBUtils::listGlobalTags()
{
  nlohmann::json resp = cdbclient->getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  for (auto &it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    std::cout << "global tag: " <<  exist_gt << std::endl;
  }
  return;
}

void CDBUtils::listPayloadTypes()
{
  nlohmann::json resp = cdbclient->getPayloadTypes();
  nlohmann::json msgcont = resp["msg"];
  for (auto &it : msgcont.items())
  {
    std::string exist_pl = it.value().at("name");
    std::cout << "payload type: " <<  exist_pl << std::endl;
  }
  return;
}

int CDBUtils::insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start)
{
  if (!isGlobalTagSet())
  {
    std::cout << "No Global Tag set" << std::endl;
    return -1;
  }
  nlohmann::json resp = cdbclient->insertPayload(pl_type,file_url,iov_start);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error inserting payload " << file_url << ", msg: " << resp["msg"] << std::endl;
  }
  else
  {
    std::cout << resp << std::endl;
  }
  return iret;
}

int CDBUtils::insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end)
{
  if (!isGlobalTagSet())
  {
    std::cout << "No Global Tag set" << std::endl;
    return -1;
  }
  nlohmann::json resp = cdbclient->insertPayload(pl_type, file_url, iov_start, iov_end);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error inserting payload " << file_url << ", msg: " << resp["msg"] << std::endl;
  }
  else
  {
    std::cout << resp << std::endl;
  }
  return iret;
}

int CDBUtils::setGlobalTag(const std::string &tagname)
{
  nlohmann::json resp = cdbclient->setGlobalTag(tagname);
  int iret = resp["code"];
  if (iret != 0)
  {
      std::cout << "Error setting global tag, msg: " << resp["msg"] << std::endl;
    }
  std::cout << "message: " << resp["msg"] << std::endl;
  return iret;
}

bool CDBUtils::isGlobalTagSet()
{
  return cdbclient->isGlobalTagSet();
}

int CDBUtils::createPayloadType(const std::string& pl_type)
{
  if (! isGlobalTagSet())
  {
    std::cout << "No Global Tag Set to add payload " << pl_type << " to" << std::endl;
    return -1;
  }
  nlohmann::json resp = cdbclient->createPayloadType(pl_type);
  nlohmann::json msgcont = resp["msg"];
  for (auto &it : msgcont.items())
  {
    std::cout << it.value() << std::endl;
    // std::string exist_pl = it.value().at("name");
    // std::cout << "payload type: " <<  exist_pl << std::endl;
  }

  int iret = 0;//resp["code"];
  if (iret != 0)
  {
      std::cout << "Error setting global tag, msg: " << resp["msg"] << std::endl;
    }
  return iret;
}
//  nlohmann::json ret = insertPayload(pl_type, file_url, 0, iov_start);
