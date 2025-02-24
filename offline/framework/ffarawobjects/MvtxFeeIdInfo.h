#ifndef FUN4ALLRAW_MVTXFEEIDINFO_H
#define FUN4ALLRAW_MVTXFEEIDINFO_H

#include <phool/PHObject.h>

#include <limits>

class MvtxFeeIdInfo : public PHObject
{
 public:
  //! ctor
  MvtxFeeIdInfo() = default;

  //! cp/mv ctor
  MvtxFeeIdInfo(const MvtxFeeIdInfo&) = default;
  MvtxFeeIdInfo(MvtxFeeIdInfo&&) = default;

  //! cp/mv assignment
  MvtxFeeIdInfo& operator=(const MvtxFeeIdInfo&) = default;
  MvtxFeeIdInfo& operator=(MvtxFeeIdInfo&&) = default;

  //! dtor
  ~MvtxFeeIdInfo() override = default;

  virtual uint16_t get_feeId() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_feeId(const uint16_t /*dummy*/) { return; }

  virtual uint32_t get_detField() const { return std::numeric_limits<uint32_t>::max(); }
  virtual void set_detField(const uint32_t /*dummy*/) { return; }

  virtual uint64_t get_bco() const { return std::numeric_limits<uint64_t>::max(); }
  virtual void set_bco(const uint64_t /*dummy*/) { return; }

 private:
  ClassDefOverride(MvtxFeeIdInfo, 1)
};

#endif  // FUN4ALLRAW_MVTXFEEIDINFO
