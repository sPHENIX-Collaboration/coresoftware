/**
 * @file trackbase/TrkrHitv3.h
 * @author Cheng-Wei Shih
 * @brief Derived class v3 for hit object with INTT timing information
 */
#ifndef TRACKBASE_TRKRHITV3_H
#define TRACKBASE_TRKRHITV3_H

#include "TrkrHitv2.h"

#include <cstdint>
#include <iostream>

class TrkrHitv3 : public TrkrHitv2
{
 public:
  //! ctor
  explicit TrkrHitv3() = default;

  //! dtor
  ~TrkrHitv3() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrHitv3 class with adc = " << m_adc
       << " and FPHX_BCO = " << m_fphx_bco
       << " and BCO = " << m_bco << std::endl;
  }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  void CopyFrom(const TrkrHit&) override;

  //! copy content from base class
  void CopyFrom(TrkrHit* source) override
  {
    CopyFrom(*source);
  }

  void setFPHXBCO(const uint16_t bco) override { m_fphx_bco = bco; }
  uint16_t getFPHXBCO() const override { return m_fphx_bco; }
  void setBCO(const uint64_t bco) override { m_bco = bco; }
  uint64_t getBCO() const override { return m_bco; }

 protected:
  uint16_t m_fphx_bco = 0;
  uint64_t m_bco = 0;

  ClassDefOverride(TrkrHitv3, 1);
};

#endif  // TRACKBASE_TRKRHITV3_H
