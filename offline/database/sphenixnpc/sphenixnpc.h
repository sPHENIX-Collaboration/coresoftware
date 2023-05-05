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
  using nopayloadclient::Client::createGlobalTag;
  using nopayloadclient::Client::deleteGlobalTag;
  using nopayloadclient::Client::setGlobalTag;

  sphenixnpc(const std::string &globaltag);
  virtual ~sphenixnpc() = default;
  nlohmann::json getUrlDict(long long iov);
  nlohmann::json getUrlDict(long long iov1, long long iov2) override;
  int createThisGlobalTag(const std::string &tagname);
  int createDomain(const std::string &domain);
  nlohmann::json get(const std::string &pl_type, long long iov);
  int cache_set_GlobalTag(const std::string &name);
  nlohmann::json clearCache() override;
  std::string getCalibrationFile(const std::string &type, uint64_t iov);
  int insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start);
  int insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end);
  int deleteThisGlobalTag(const std::string &);

  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  nlohmann::json url_dict_;  // valid until global tag is switched
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
};

#endif  // SPHENIXNPC_SPHENIXNPC_H
