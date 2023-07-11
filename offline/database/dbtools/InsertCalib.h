#ifndef DBTOOLS_INSERTCALIB_H
#define DBTOOLS_INSERTCALIB_H

#include <cstdint>
#include <string>

class InsertCalib
{
 public:
  InsertCalib() = default;
  ~InsertCalib() = default;
  int insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start);
  int insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start, uint64_t iov_end);
  int createGlobalTag(const std::string &tag);
  int listGlobalTags();
  int deleteGlobalTag(const std::string &tag);
  int createPayloadType(const std::string &payloadtype);
  int listPayloadTypes();
};

#endif
