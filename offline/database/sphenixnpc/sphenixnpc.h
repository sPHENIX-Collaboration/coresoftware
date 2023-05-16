#ifndef SPHENIXNPC_SPHENIXNPC_H
#define SPHENIXNPC_SPHENIXNPC_H

#include <nopayloadclient/nopayloadclient.hpp>

#include <nlohmann/json.hpp>

#include <cstdint>
#include <iostream>
#include <set>
#include <string>

class sphenixnpc : public nopayloadclient::Client
{
 public:
  sphenixnpc();
  sphenixnpc(const std::string &globaltag);
  virtual ~sphenixnpc() = default;
  nlohmann::json getUrlDict(long long iov);
  nlohmann::json getUrlDict(long long iov1, long long iov2) override;
  nlohmann::json createGlobalTag(const std::string &tagname) override;
  int createDomain(const std::string &domain);
  nlohmann::json setGlobalTag(const std::string &tagname) override;
  nlohmann::json get(const std::string &pl_type, long long iov);
  int cache_set_GlobalTag(const std::string &name);
  nlohmann::json clearCache() override;
  std::string getUrl(const std::string &type, uint64_t iov);
  nlohmann::json insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start);
  nlohmann::json insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end);
  bool isGlobalTagSet();
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  uint64_t m_CachedIOV = 0;
  nlohmann::json url_dict_;  // valid until global tag is switched
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
};

#endif  // SPHENIXNPC_SPHENIXNPC_H
