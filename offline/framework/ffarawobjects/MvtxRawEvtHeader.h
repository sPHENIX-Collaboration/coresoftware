#ifndef FUN4ALLRAW_MVTXRAWEVTHEADER_H
#define FUN4ALLRAW_MVTXRAWEVTHEADER_H

#include <phool/PHObject.h>

#include <cstdint>
#include <limits>
#include <map>
#include <set>

class MvtxFeeIdInfo;

class MvtxRawEvtHeader : public PHObject
{
 public:
  //! ctor
  MvtxRawEvtHeader() = default;

  //! cp/mv ctor
  MvtxRawEvtHeader(const MvtxRawEvtHeader &) = default;
  MvtxRawEvtHeader(MvtxRawEvtHeader &&) = default;

  //! cp/mv assignment
  MvtxRawEvtHeader &operator=(const MvtxRawEvtHeader &) = default;
  MvtxRawEvtHeader &operator=(MvtxRawEvtHeader &&) = default;

  //! dtor
  ~MvtxRawEvtHeader() override = default;

  virtual void AddFeeId(const int & /*dummy*/) { return; }
  virtual void AddL1Trg(const uint64_t & /*dummy*/) { return; }

  virtual void AddFeeId(const std::set<uint16_t> & /*dummy*/) { return; }
  virtual void AddL1Trg(const std::set<uint64_t> & /*dummy*/) { return; }

  virtual std::set<uint16_t> &getMvtxFeeIdSet() { return dummySet16; }
  virtual std::set<uint64_t> &getMvtxLvL1BCO() { return dummySet64; }

  virtual MvtxFeeIdInfo *AddFeeIdInfo() { return nullptr; }
  virtual MvtxFeeIdInfo *AddFeeIdInfo(MvtxFeeIdInfo * /*dummy*/) { return nullptr; }

  virtual uint64_t get_nFeeIdInfo() { return 100; }
  virtual MvtxFeeIdInfo *get_feeIdInfo(unsigned int /*dummy*/) { return nullptr; }

 private:
  std::set<uint16_t> dummySet16;
  std::set<uint64_t> dummySet64;

  ClassDefOverride(MvtxRawEvtHeader, 1)
};

#endif
