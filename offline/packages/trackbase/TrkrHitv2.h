/**
 * @file trackbase/TrkrHit.h
 * @author D. McGlinchey
 * @date 4 June 2018
 * @brief Base class for hit object
 */
#ifndef TRACKBASE_TRKRHITV2_H
#define TRACKBASE_TRKRHITV2_H

#include "TrkrDefs.h"
#include "TrkrHit.h"

#include <phool/PHObject.h>

#include <iostream>

/**
 * @brief Base class for hit object
 *
 * This is the empyt virtual base class for a hit object.
 * Each subsystem should implement an inherited version
 * which contains the actual storage information.
 */
class TrkrHitv2 : public TrkrHit
{
 public:
  //! ctor
  TrkrHitv2(); 

  //! dtor
  ~TrkrHitv2() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrHitv2 class with adc = " << m_adc << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  // these set and get the energy before digitization
  void addEnergy(const double edep) override;
  double getEnergy() override;

  // after digitization, these are the adc values
  void setAdc(const unsigned int adc) override;
  unsigned int getAdc() override ;

 protected:

  unsigned short m_adc = 0;
  ClassDefOverride(TrkrHitv2, 1);
};

#endif //TRACKBASE_TRKRHITV2_H
