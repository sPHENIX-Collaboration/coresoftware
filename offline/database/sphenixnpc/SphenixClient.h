#ifndef SPHENIXNPC_SPHENIXCLIENT_H
#define SPHENIXNPC_SPHENIXCLIENT_H

#include <nopayloadclient/nopayloadclient.hpp>

#include <nlohmann/json.hpp>

#include <cstdint>
#include <iostream>
#include <set>
#include <string>

class SphenixClient : public nopayloadclient::NoPayloadClient
{
 public:
  SphenixClient() = default;
  SphenixClient(const std::string& globaltag);
  virtual ~SphenixClient() = default;
  // make clang happy, since we use our own without overriding the base class methods
  using nopayloadclient::NoPayloadClient::getPayloadIOVs;

  nlohmann::json getPayloadIOVs(long long iov);
  nlohmann::json getUrl(const std::string& pl_type, long long iov);
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start);
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start, long long iov_end);
  nlohmann::json setGlobalTag(const std::string& name);
  nlohmann::json clearCache();
  std::string getCalibration(const std::string& pl_type, long long iov);

  int createDomain(const std::string& domain);
  int cache_set_GlobalTag(const std::string& name);
  bool isGlobalTagSet();
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
};

#endif  // SPHENIXNPC_SPHENIXCLIENT_H
