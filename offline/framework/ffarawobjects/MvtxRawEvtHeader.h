#ifndef FUN4ALLRAW_MVTXRAWEVTHEADER_H
#define FUN4ALLRAW_MVTXRAWEVTHEADER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <limits>
#include <set>
#include <map>

class MvtxFeeIdInfo;

class MvtxRawEvtHeader : public PHObject
{
 public:
  MvtxRawEvtHeader() = default;
  virtual ~MvtxRawEvtHeader() = default;

  virtual void AddFeeId(const int&) { return; }
  virtual void AddL1Trg(const uint64_t&) { return; }

  virtual void AddFeeId(const std::set<uint16_t>&) { return; }
  virtual void AddL1Trg(const std::set<uint64_t>&) { return; }

  virtual std::set<uint16_t>& getMvtxFeeIdSet() { return dummySet16; }
  virtual std::set<uint64_t>& getMvtxLvL1BCO() { return dummySet64; }

  virtual MvtxFeeIdInfo *AddFeeIdInfo() { return nullptr; }
  virtual MvtxFeeIdInfo *AddFeeIdInfo(MvtxFeeIdInfo *) { return nullptr; }

  virtual uint64_t get_nFeeIdInfo() { return 100; }
  virtual MvtxFeeIdInfo *get_feeIdInfo(unsigned int) { return nullptr; }


 private:
  std::set<uint16_t> dummySet16;
  std::set<uint64_t> dummySet64;

  ClassDefOverride(MvtxRawEvtHeader, 1)
};

#endif
