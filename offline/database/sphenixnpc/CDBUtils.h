#ifndef SPHENIXNPC_CDBUTILS_H
#define SPHENIXNPC_CDBUTILS_H

#include <cstdint>
#include <iostream>
#include <set>
#include <string>

class sphenixnpc;

class CDBUtils
{
 public:

  CDBUtils();
  CDBUtils(const std::string &globaltag);
  virtual ~CDBUtils() = default;
  void createGlobalTag(const std::string &tagname);
  void setGlobalTag(const std::string &tagname);
  int createDomain(const std::string &domain);
  std::string getCalibrationFile(const std::string &type, uint64_t iov);
  int insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start);
  int insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end);
  int deleteGlobalTag(const std::string &);
  void ListGlobalTags();
  void clearCache();
  void Verbosity(int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 private:
  int m_Verbosity = 0;
  sphenixnpc *cdbclient = nullptr;
  std::string m_CachedGlobalTag;
  std::set<std::string> m_DomainCache;
};

#endif  // SPHENIXNPC_SPHENIXNPC_H
