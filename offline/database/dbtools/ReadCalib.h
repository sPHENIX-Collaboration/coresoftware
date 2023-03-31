#ifndef DBTOOLS_READCALIB_H
#define DBTOOLS_READCALIB_H

#include <cstdint>
#include <string>

class ReadCalib
{
 public:
  ReadCalib() = default;
  ~ReadCalib() = default;
  int listGlobalTags();
  int listPayloadTypes();
  std::string getCalibrationFile(const std::string &type, uint64_t iov) const;
  std::string getCalibrationFile(const std::string &globaltag, const std::string &type, uint64_t iov) const;
};

#endif
