#ifndef SPHENIXNPC_CDBUTILS_H
#define SPHENIXNPC_CDBUTILS_H

#include <cstdint>  // for uint64_t
#include <set>
#include <string>

class SphenixClient;

class CDBUtils
{
 public:
  CDBUtils();
  explicit CDBUtils(const std::string &globaltag);

  // delete copy ctor and assignment operator (cppcheck)
  explicit CDBUtils(const CDBUtils &) = delete;
  CDBUtils &operator=(const CDBUtils &) = delete;

  virtual ~CDBUtils() = default;
  int createGlobalTag(const std::string &tagname);
  int setGlobalTag(const std::string &tagname);
  int lockGlobalTag(const std::string &tagname);
  int unlockGlobalTag(const std::string &tagname);
  int createPayloadType(const std::string &domain);
  std::string getUrl(const std::string &type, uint64_t iov);
  int insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start);
  int insertPayload(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end);
  int cloneGlobalTag(const std::string &source, const std::string &target);

  int deleteGlobalTag(const std::string &);
  void listGlobalTags();
  void listPayloadTypes();
  void listPayloadIOVs(uint64_t iov);
  void clearCache();
  bool isGlobalTagSet();
  void Verbosity(int i);
  int Verbosity() const { return m_Verbosity; }
  int deletePayloadIOV(const std::string &pl_type, uint64_t iov_start);
  int deletePayloadIOV(const std::string &pl_type, uint64_t iov_start, uint64_t iov_end);

 private:
  int m_Verbosity = 0;
  SphenixClient *cdbclient = nullptr;
  std::string m_CachedGlobalTag;
  std::set<std::string> m_PayloadTypeCache;
};

#endif  // SPHENIXNPC_CDBUTILS_H
