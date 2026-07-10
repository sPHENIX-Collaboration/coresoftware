/**
 * @file trackbase/TrkrHitv1.h
 * @author D. McGlinchey, Tony Frawley
 * @brief Derived class v1 for hit object
 */
#ifndef TRACKBASE_TRKRHITV1_H
#define TRACKBASE_TRKRHITV1_H

#include "TrkrHit.h"

#include <phool/PHObject.h>

#include <iostream>

/**
 * @brief Inherited class v1 for hit object
 *
 */
class TrkrHitv1 : public TrkrHit
{
 public:
  //! ctor
  explicit TrkrHitv1() = default;

  //! dtor
  ~TrkrHitv1() override = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrHitV1 class with adc = " << m_adc << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  void CopyFrom(const TrkrHit&) override;

  //! copy content from base class
  void CopyFrom(TrkrHit* source) override
  {
    CopyFrom(*source);
  }

  void addEnergy(const double edep) override { m_edep += edep; }
  double getEnergy() const override { return m_edep; }
  void setAdc(const unsigned int adc) override { m_adc = adc; }
  unsigned int getAdc() const override;

 protected:
  double m_edep = 0;
  unsigned int m_adc = 0;
  ClassDefOverride(TrkrHitv1, 1);
};

#endif  // TRACKBASE_TRKRHITV1_H
