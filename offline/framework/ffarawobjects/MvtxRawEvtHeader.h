#ifndef FUN4ALLRAW_MVTXRAWEVTHEADER_H
#define FUN4ALLRAW_MVTXRAWEVTHEADER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <limits>
#include <set>
#include <map>

class MvtxRawEvtHeader : public PHObject
{
 public:
  MvtxRawEvtHeader() = default;
  virtual ~MvtxRawEvtHeader() = default;

  virtual void AddFeeId(const int&) { return; }
  virtual void AddL1Trg(const uint64_t&) { return; }
  virtual void AddDetField(const uint16_t&,
                           const uint32_t&) { return; }

  virtual void AddFeeId(const std::set<uint16_t>&) { return; }
  virtual void AddL1Trg(const std::set<uint64_t>&) { return; }
  virtual void AddDetField(const std::map<uint16_t, uint32_t>&) { return; }

  virtual std::set<uint16_t>& getMvtxFeeIdSet() { return dummySet16; }
  virtual std::set<uint64_t>& getMvtxLvL1BCO() { return dummySet64; }
  virtual std::map<uint16_t, uint32_t>& getMvtxDetField() { return dummyMap; }

 private:
  std::set<uint16_t> dummySet16;
  std::set<uint64_t> dummySet64;
  std::map<uint16_t, uint32_t> dummyMap;

  ClassDefOverride(MvtxRawEvtHeader, 2)
};

#endif
