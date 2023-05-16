#include "sphenixnpc.h"

#include <nopayloadclient/nopayloadclient.hpp>

#include <nlohmann/json.hpp>
#include <iostream>
#include <stdexcept>

sphenixnpc::sphenixnpc()
{}

sphenixnpc::sphenixnpc(const std::string &globaltag)
{
  cache_set_GlobalTag(globaltag);
}


nlohmann::json sphenixnpc::createGlobalTag(const std::string &tagname)
{
  if (tagname == m_CachedGlobalTag) // global tag already set
  {
    std::string message = "global tag already set to " + tagname;
    return {{"code", 0}, {"msg", message}};
  }
// check if the global tag exists already
  bool found_gt = false;
  nlohmann::json resp = nopayloadclient::Client::getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  for (auto &it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    std::cout << "global tag: " <<  exist_gt << std::endl;
    if (exist_gt == tagname)
    {
      found_gt = true;
      break;
    }
  }
  if (found_gt)
  {
    std::string message = "global tag " + tagname + " already exists";
    return {{"code", 0}, {"msg", message}};
  }
  resp = nopayloadclient::Client::createGlobalTag(tagname);
  int iret = resp["code"];
  if (iret != 0)
  {
    std::cout << "Error creating global tag, msg: " << resp["msg"] << std::endl;
  }
  else
  {
    nopayloadclient::Client::setGlobalTag(tagname);
    m_CachedGlobalTag = tagname;
    clearCache();
    std::cout << "sphenixnpc: Created new global tag " << tagname << std::endl;
  }
  return resp;
}

nlohmann::json sphenixnpc::getUrlDict(long long iov)
{
  return nopayloadclient::Client::getUrlDict(iov,iov);
}

nlohmann::json sphenixnpc::getUrlDict(long long iov1, long long iov2)
{
  return nopayloadclient::Client::getUrlDict(iov1,iov2);
}

nlohmann::json sphenixnpc::get(const std::string &pl_type, long long iov)
{
  if (url_dict_.is_null())
  {
    nlohmann::json resp = getUrlDict(iov);
    std::cout << "response" << std::endl;
    std::cout << resp << std::endl;
    std::cout << "end response" << std::endl;
    if (resp["code"] != 0)
    {
      return resp;
    }
    url_dict_ = resp["msg"];
  }
  if (not url_dict_.contains(pl_type))
  {
    return nopayloadclient::DataBaseException("No payload with type " + pl_type + " exists.").jsonify();
  }
  return makeResp(url_dict_[pl_type]);
}

nlohmann::json sphenixnpc::setGlobalTag(const std::string &tagname)
{
  if (tagname == m_CachedGlobalTag) // global tag already set
  {
    std::string message = "global tag already set to " + tagname;
    return {{"code", 0}, {"msg", message}};
  }
  url_dict_ = nlohmann::json{};
  bool found_gt = false;
  nlohmann::json resp = nopayloadclient::Client::getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  for (auto &it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    std::cout << "global tag: " <<  exist_gt << std::endl;
    if (exist_gt == tagname)
    {
      found_gt = true;
      break;
    }
  }
  if (found_gt)
  {
  m_CachedGlobalTag = tagname;
  resp = nopayloadclient::Client::setGlobalTag(tagname);
  return resp;
  }
std::string message = "global tag " + tagname + " does not exist";
    return {{"code", -1}, {"msg", message}};
}

int sphenixnpc::cache_set_GlobalTag(const std::string &tagname)
{
  int iret = 0;
  if (tagname == m_CachedGlobalTag) // global tag already set
  {
    return iret;
  }
  url_dict_ = nlohmann::json{};
  m_CachedGlobalTag = tagname;
  nopayloadclient::Client::setGlobalTag(tagname);
  bool found_gt = false;
  nlohmann::json resp = nopayloadclient::Client::getGlobalTags();
  nlohmann::json msgcont = resp["msg"];
  for (auto &it : msgcont.items())
  {
    std::string exist_gt = it.value().at("name");
    std::cout << "global tag: " <<  exist_gt << std::endl;
    if (exist_gt == tagname)
    {
      found_gt = true;
      break;
    }
  }
  if (!found_gt)
  {
    resp = nopayloadclient::Client::createGlobalTag();
    iret = resp["code"];
    if (iret != 0)
    {
      std::cout << "Error creating global tag, msg: " << resp["msg"] << std::endl;
    }
    else
    {
      std::cout << "sphenixnpc: Created new global tag " << tagname << std::endl;
    }
  }
  return iret;
}

nlohmann::json sphenixnpc::clearCache()
{
  url_dict_ = nlohmann::json{};
  return nopayloadclient::Client::clearCache();
}

std::string sphenixnpc::getUrl(const std::string &type, uint64_t iov)
{
  if (iov != m_CachedIOV)
  {
    clearCache();
    m_CachedIOV = iov;
  }
  nlohmann::json result = get(type, iov);
  std::string fullUrl = result.at("msg");
  if (fullUrl.find("Exception") != std::string::npos)
  {
    if (Verbosity())
    {
      std::cout << fullUrl << std::endl;
    }
    return std::string("");
  }
  return fullUrl;
}

nlohmann::json sphenixnpc::insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start)
{
  nlohmann::json resp = nopayloadclient::Client::insertPayload(pl_type, file_url, 0, iov_start);
  if (Verbosity())
  {
    std::cout << resp << std::endl;
  }
  return resp;
}

nlohmann::json sphenixnpc::insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end)
{
  nlohmann::json ret = nopayloadclient::Client::insertPayload(pl_type, file_url, 0, iov_start, 0, iov_end);
  if (Verbosity())
  {
    std::cout << ret << std::endl;
  }
  return ret;
}

int sphenixnpc::createDomain(const std::string &domain)
{
  int iret = -1;
  nlohmann::json resp;
  if (m_DomainCache.empty())
  {
    resp = nopayloadclient::Client::getPayloadTypes();
    nlohmann::json msgcont = resp["msg"];
    for (auto &it : msgcont.items())
    {
      std::string existent_domain = it.value().at("name");
      m_DomainCache.insert(existent_domain);
    }
  }
  if (m_DomainCache.find(domain) == m_DomainCache.end())
  {
    resp = nopayloadclient::Client::createPayloadType(domain);
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

bool sphenixnpc::isGlobalTagSet()
{
  if (m_CachedGlobalTag.empty())
  {
    return false;
  }
  return true;
}

