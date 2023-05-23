#ifndef SPHENIXNPC_CDBUTILS_H
#define SPHENIXNPC_CDBUTILS_H

#include <cstdint>
#include <iostream>
#include <set>
#include <string>

class SphenixClient;

class CDBUtils
{
 public:

  CDBUtils();
  CDBUtils(const std::string &globaltag);
  virtual ~CDBUtils() = default;
  void createGlobalTag(const std::string &tagname);
  int setGlobalTag(const std::string &tagname);
  void lockGlobalTag(const std::string &tagname);
  void unlockGlobalTag(const std::string &tagname);
  int createDomain(const std::string &domain);
  std::string getUrl(const std::string &type, uint64_t iov);
  int insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start);
  int insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end);
  int deleteGlobalTag(const std::string &);
  void listGlobalTags();
  void listPayloadTypes();
  void clearCache();
  bool isGlobalTagSet();
  int createPayloadType(const std::string& pl_type);
  void Verbosity(int i);
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  SphenixClient *cdbclient = nullptr;
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
};

#endif  // SPHENIXNPC_CDBUTILS_H
