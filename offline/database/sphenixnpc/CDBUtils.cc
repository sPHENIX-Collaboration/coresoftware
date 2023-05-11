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
  cdbclient->createThisGlobalTag(tagname);
}

void CDBUtils::setGlobalTag(const std::string &tagname)
{
  cdbclient->setGlobalTag(const std::string &tagname);
}

int CDBUtils::deleteGlobalTag(const std::string &tagname)
{
  int iret = cdbclient->deleteThisGlobalTag(tagname);
  return iret;
}


void CDBUtils::clearCache()
{
  cdbclient->clearCache();
}

std::string CDBUtils::getCalibrationFile(const std::string &type, uint64_t iov)
{
  return cdbclient->getCalibrationFile(type, iov);
}

int CDBUtils::insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start)
{
  int iret = cdbclient->insertcalib(pl_type, file_url, iov_start);
  return iret;
}

int CDBUtils::insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end)
{
  int iret = cdbclient->insertcalib(pl_type, file_url, iov_start,iov_end);
  return iret;
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

void CDBUtils::ListGlobalTags()
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

