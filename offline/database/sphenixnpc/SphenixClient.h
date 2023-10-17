#ifndef SPHENIXNPC_SPHENIXCLIENT_H
#define SPHENIXNPC_SPHENIXCLIENT_H

#include <nopayloadclient/nopayloadclient.hpp>

#include <nlohmann/json.hpp>

#include <set>
#include <string>

class SphenixClient : public nopayloadclient::NoPayloadClient
{
 public:
  SphenixClient() = default;
  explicit SphenixClient(const std::string& globaltag);
  virtual ~SphenixClient() = default;
  // make clang happy, since we use our own without overriding the base class methods
  using nopayloadclient::NoPayloadClient::getPayloadIOVs;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start);
  nlohmann::json getUrlDict(long long iov);
  nlohmann::json deletePayloadIOV(const std::string& pl_type, long long iov_start);
#pragma GCC diagnostic pop

  nlohmann::json getPayloadIOVs(long long iov);
  nlohmann::json getUrl(const std::string& pl_type, long long iov);
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start, long long iov_end) override;
  nlohmann::json setGlobalTag(const std::string& name) override;
  std::string getCalibration(const std::string& pl_type, long long iov);
  nlohmann::json unlockGlobalTag(const std::string& tagname) override;
  nlohmann::json lockGlobalTag(const std::string& tagname) override;
  nlohmann::json deletePayloadIOV(const std::string& pl_type, long long iov_start, long long iov_end) override;

  bool existGlobalTag(const std::string& tagname);
  int createDomain(const std::string& domain);
  int cache_set_GlobalTag(const std::string& name);
  bool isGlobalTagSet();
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
  std::set<std::string> m_GlobalTagCache;
};

#endif  // SPHENIXNPC_SPHENIXCLIENT_H
