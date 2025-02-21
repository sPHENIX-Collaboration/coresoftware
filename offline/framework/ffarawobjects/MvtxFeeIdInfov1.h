#ifndef FUN4ALLRAW_MVTXFEEIDINFOV1_H
#define FUN4ALLRAW_MVTXFEEIDINFOV1_H

#include "MvtxFeeIdInfo.h"

#include <limits>

class MvtxFeeIdInfov1 : public MvtxFeeIdInfo
{
 public:
  //! ctor
  MvtxFeeIdInfov1() = default;
  explicit MvtxFeeIdInfov1(MvtxFeeIdInfo* info);

  //! cp/mv ctor
  MvtxFeeIdInfov1(const MvtxFeeIdInfov1&) = default;
  MvtxFeeIdInfov1(MvtxFeeIdInfov1&&) = default;

  //! cp/mv assigment
  MvtxFeeIdInfov1& operator=(const MvtxFeeIdInfov1&) = default;
  MvtxFeeIdInfov1& operator=(MvtxFeeIdInfov1&&) = default;

  //! dtor
  ~MvtxFeeIdInfov1() override = default;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  uint16_t get_feeId() const override { return m_feeId; }
  // cppcheck-suppress virtualCallInConstructor
  void set_feeId(const uint16_t val) override { m_feeId = val; }

  uint32_t get_detField() const override { return m_detField; }
  // cppcheck-suppress virtualCallInConstructor
  void set_detField(const uint32_t val) override { m_detField = val; }

  uint64_t get_bco() const override { return m_bco; }
  // cppcheck-suppress virtualCallInConstructor
  void set_bco(const uint64_t val) override { m_bco = val; }

 private:
  uint16_t m_feeId = std::numeric_limits<uint16_t>::max();
  uint32_t m_detField = std::numeric_limits<uint32_t>::max();
  uint64_t m_bco = std::numeric_limits<uint64_t>::max();

  ClassDefOverride(MvtxFeeIdInfov1, 1)
};

#endif
