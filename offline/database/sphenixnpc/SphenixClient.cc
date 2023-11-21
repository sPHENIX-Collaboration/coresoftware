#include "SphenixClient.h"

#include <nopayloadclient/exception.hpp>  // for DataBaseException
#include <nopayloadclient/nopayloadclient.hpp>

#include <nlohmann/json.hpp>

#include <iostream>
#include <stdexcept>

SphenixClient::SphenixClient(const std::string& gt_name)
  : nopayloadclient::NoPayloadClient(gt_name)
  , m_CachedGlobalTag(gt_name)
{
}

nlohmann::json SphenixClient::getPayloadIOVs(long long iov)
{
  return nopayloadclient::NoPayloadClient::getPayloadIOVs(0, iov);
}

nlohmann::json SphenixClient::getUrl(const std::string& pl_type, long long iov)
{
  nlohmann::json resp = getPayloadIOVs(iov);
  if (resp["code"] != 0)
  {
    return resp;
  }
  nlohmann::json payload_iovs = resp["msg"];
  if (!payload_iovs.contains(pl_type))
  {
    return nopayloadclient::DataBaseException("No valid payload with type " + pl_type).jsonify();
  }
  nlohmann::json payload_iov = payload_iovs[pl_type];
  if (m_Verbosity > 0)
  {
    std::cout << "pl_type: " << pl_type
              << ", iov: " << iov
              << ", minor_iov_start: " << payload_iov["minor_iov_start"]
              << ", minor_iov_end: " << payload_iov["minor_iov_end"]
              << std::endl;
  }
  if (payload_iov["minor_iov_end"] <= iov)
  {
    return nopayloadclient::DataBaseException("No valid payload with type " + pl_type).jsonify();
  }
  return makeResp(payload_iov["payload_url"]);
}

nlohmann::json SphenixClient::getUrlDict(long long iov)
{
  nlohmann::json resp = getPayloadIOVs(iov);
  if (resp["code"] != 0)
  {
    return resp;
  }
  for (auto it = resp["msg"].begin(); it != resp["msg"].end();)
  {
    if (it.value()["minor_iov_end"] < iov)
    {
      it = resp["msg"].erase(it);
    }
    else
    {
      ++it;
    }
  }
  for (auto& piov : resp["msg"].items())
  {
    piov.value() = piov.value()["payload_url"];
  }
  return resp;
}

nlohmann::json SphenixClient::deletePayloadIOV(const std::string& pl_type, long long iov_start)
{
  return nopayloadclient::NoPayloadClient::deletePayloadIOV(pl_type, 0, iov_start);
}

nlohmann::json SphenixClient::deletePayloadIOV(const std::string& pl_type, long long iov_start, long long iov_end)
{
  return nopayloadclient::NoPayloadClient::deletePayloadIOV(pl_type, 0, iov_start, 0, iov_end);
}

std::string SphenixClient::getCalibration(const std::string& pl_type, long long iov)
{
  nlohmann::json resp = getUrl(pl_type, iov);
  if (resp["code"] != 0)
  {
    if (m_Verbosity > 0)
    {
      std::cout << resp << std::endl;
    }
    return "";
  }
  return resp["msg"];
}

nlohmann::json SphenixClient::unlockGlobalTag(const std::string& gt_name)
{
  if (existGlobalTag(gt_name))
  {
    return nopayloadclient::NoPayloadClient::unlockGlobalTag(gt_name);
  }
  std::string message = "global tag " + gt_name + " does not exist";
  return {{"code", -1}, {"msg", message}};
}

nlohmann::json SphenixClient::lockGlobalTag(const std::string& gt_name)
{
  if (existGlobalTag(gt_name))
  {
    return nopayloadclient::NoPayloadClient::lockGlobalTag(gt_name);
  }
  std::string message = "global tag " + gt_name + " does not exist";
  return {{"code", -1}, {"msg", message}};
}

nlohmann::json SphenixClient::insertPayload(const std::string& pl_type, const std::string& file_url,
                                            long long iov_start)
{
  return nopayloadclient::NoPayloadClient::insertPayload(pl_type, file_url, 0, iov_start);
}

nlohmann::json SphenixClient::insertPayload(const std::string& pl_type, const std::string& file_url,
                                            long long iov_start, long long iov_end)
{
  return nopayloadclient::NoPayloadClient::insertPayload(pl_type, file_url, 0, iov_start, 0, iov_end);
}

nlohmann::json SphenixClient::setGlobalTag(const std::string& gt_name)
{
  if (m_CachedGlobalTag == gt_name)
  {
    std::string message = "global tag already set to " + gt_name;
    return {{"code", 0}, {"msg", message}};
  }
  // check if the global tag exists before switching
  if (existGlobalTag(gt_name))
  {
    m_CachedGlobalTag = gt_name;
    return nopayloadclient::NoPayloadClient::setGlobalTag(gt_name);
  }

  std::string message = "global tag " + gt_name + " does not exist";
  return {{"code", -1}, {"msg", message}};
}

int SphenixClient::cache_set_GlobalTag(const std::string& tagname)
{
  int iret = 0;
  if (tagname == m_CachedGlobalTag)  // global tag already set
  {
    return iret;
  }
  m_CachedGlobalTag = tagname;
  nopayloadclient::NoPayloadClient::setGlobalTag(tagname);
  bool found_gt = false;
  nlohmann::json resp = nopayloadclient::NoPayloadClient::getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  for (auto& it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    std::cout << "global tag: " << exist_gt << std::endl;
    if (exist_gt == tagname)
    {
      found_gt = true;
      break;
    }
  }
  if (!found_gt)
  {
    resp = nopayloadclient::NoPayloadClient::createGlobalTag();
    iret = resp["code"];
    if (iret != 0)
    {
      std::cout << "Error creating global tag, msg: " << resp["msg"] << std::endl;
    }
    else
    {
      std::cout << "SphenixClient: Created new global tag " << tagname << std::endl;
    }
  }
  return iret;
}

int SphenixClient::createDomain(const std::string& domain)
{
  int iret = -1;
  nlohmann::json resp;
  if (m_DomainCache.empty())
  {
    resp = nopayloadclient::NoPayloadClient::getPayloadTypes();
    nlohmann::json msgcont = resp["msg"];
    for (auto& it : msgcont.items())
    {
      std::string existent_domain = it.value().at("name");
      m_DomainCache.insert(existent_domain);
    }
  }
  if (m_DomainCache.find(domain) == m_DomainCache.end())
  {
    resp = nopayloadclient::NoPayloadClient::createPayloadType(domain);
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

bool SphenixClient::isGlobalTagSet()
{
  if (m_CachedGlobalTag.empty())
  {
    return false;
  }
  return true;
}

bool SphenixClient::existGlobalTag(const std::string& tagname)
{
  if (m_GlobalTagCache.find(tagname) != m_GlobalTagCache.end())
  {
    return true;
  }
  nlohmann::json resp = nopayloadclient::NoPayloadClient::getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  for (auto& it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    m_GlobalTagCache.insert(tagname);
    if (exist_gt == tagname)
    {
      return true;
    }
  }
  return false;
}
