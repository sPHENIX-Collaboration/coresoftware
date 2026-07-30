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
  explicit SphenixClient(const std::string& gt_name);
  virtual ~SphenixClient() = default;
  // Keep the base overloads visible while adding sPHENIX convenience wrappers.
  using nopayloadclient::NoPayloadClient::deletePayloadIOV;
  using nopayloadclient::NoPayloadClient::getPayloadIOVs;
  using nopayloadclient::NoPayloadClient::lockGlobalTag;
  using nopayloadclient::NoPayloadClient::unlockGlobalTag;

  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start);
  nlohmann::json getUrlDict(long long iov);
  nlohmann::json getUrlDict(long long major_iov, long long minor_iov) override;
  nlohmann::json deletePayloadIOV(const std::string& pl_type, long long iov_start);

  nlohmann::json getPayloadIOVs(long long iov);
  nlohmann::json getUrl(const std::string& pl_type, long long iov);
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long iov_start, long long iov_end) override;
  nlohmann::json insertPayload(const std::string& pl_type, const std::string& file_url, long long major_iov_start, long long minor_iov_start, long long major_iov_end, long long minor_iov_end) override;
  nlohmann::json setGlobalTag(const std::string& name) override;
  std::string getCalibration(const std::string& pl_type, long long iov);
  nlohmann::json unlockGlobalTag(const std::string& gt_name) override;
  nlohmann::json lockGlobalTag(const std::string& gt_name) override;
  nlohmann::json deletePayloadIOV(const std::string& pl_type, long long iov_start, long long iov_end) override;

  bool existGlobalTag(const std::string& gt_name);
  int createDomain(const std::string& domain);
  int cache_set_GlobalTag(const std::string& name);
  bool isGlobalTagSet();
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }
  void DumpCalibrations(long long iov, const std::string& filename);

 private:
  int m_Verbosity{0};
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
  std::set<std::string> m_GlobalTagCache;
};

#endif  // SPHENIXNPC_SPHENIXCLIENT_H
