/**
 * @file trackbase/RawHit.h
 * @author D. McGlinchey
 * @date 4 June 2018
 * @brief Base class for hit object
 */
#ifndef TRACKBASE_RAWHITV1_H
#define TRACKBASE_RAWHITV1_H

#include "TrkrDefs.h"
#include "RawHit.h"

#include <phool/PHObject.h>

#include <iostream>

/**
 * @brief Base class for hit object
 *
 * This is the empyt virtual base class for a hit object.
 * Each subsystem should implement an inherited version
 * which contains the actual storage information.
 */
class RawHitv1 : public RawHit
{
 public:
  //! ctor
  RawHitv1(); 

  //! dtor
  ~RawHitv1() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "RawHitv1 class with adc = " << m_adc << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  // after digitization, these are the adc values
  void setAdc(const unsigned int adc) override;
  unsigned int getAdc() override ;

  void setPhiBin(const unsigned int phibin) override;
  unsigned int getPhiBin() override;

  void setTBin(const unsigned int tbin ) override;
  unsigned int getTBin() override;

  unsigned short m_adc = 0;
  unsigned short m_tbin = 0;
  unsigned short m_phibin = 0;

 protected:

  ClassDefOverride(RawHitv1, 1);
};

#endif //TRACKBASE_RAWHITV1_H
