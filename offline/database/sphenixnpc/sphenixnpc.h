#ifndef SPHENIXNPC_SPHENIXNPC_H
#define SPHENIXNPC_SPHENIXNPC_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <nopayloadclient/nopayloadclient.hpp>
#pragma GCC diagnostic pop

#include <nlohmann/json.hpp>

#include <iostream>

class sphenixnpc : public nopayloadclient::Client
{
 public:
  static sphenixnpc *instance(const std::string &globaltag = "NONE");
  ~sphenixnpc();
  nlohmann::json getUrlDict(long long iov);
  nlohmann::json createGlobalTag(const std::string &tagname);

  nlohmann::json deleteGlobalTag(const std::string &);
  nlohmann::json get(std::string pl_type, long long iov);
  nlohmann::json insertPayload(std::string pl_type, std::string file_url,
                               long long iov_start);
  nlohmann::json insertPayload(std::string pl_type, std::string file_url,
                               long long iov_start, long long iov_end) override;
  nlohmann::json setGlobalTag(std::string name) override;
  nlohmann::json clearCache() override;
  std::string getCalibrationFile(const std::string &type, uint64_t iov);
  int insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start);
  int insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start, uint64_t iov_end);
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  static sphenixnpc *__instance;
  int m_Verbosity = 0;
  nlohmann::json url_dict_;  // valid until global tag is switched
  std::string m_CachedGlobalTag;
};

#endif  // SPHENIXNPC_SPHENIXNPC_H
