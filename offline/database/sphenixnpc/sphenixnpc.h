#ifndef SPHENIXNPC_SPHENIXNPC_H
#define SPHENIXNPC_SPHENIXNPC_H

#include <nopayloadclient/nopayloadclient.hpp>

#include <nlohmann/json.hpp>


#include <cstdint>
#include <iostream>
#include <set>
#include <string>

class sphenixnpc : public nopayloadclient::NoPayloadClient
{
 public:
  sphenixnpc();
  sphenixnpc(const std::string &globaltag);
  virtual ~sphenixnpc() = default;
  nlohmann::json getPayloadIOVs(long long iov);
  nlohmann::json  getUrl(const std::string& pl_type, long long iov);
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start);
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start, long long iov_end);
  nlohmann::json setGlobalTag1(const std::string& name);
  nlohmann::json clearCache();

  nlohmann::json createGlobalTag1(const std::string &tagname);
  int createDomain(const std::string &domain);
  nlohmann::json setGlobalTag(const std::string &tagname);
  int cache_set_GlobalTag(const std::string &name);
  bool isGlobalTagSet();
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  uint64_t m_CachedIOV = 0;
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
};

#endif  // SPHENIXNPC_SPHENIXNPC_H
