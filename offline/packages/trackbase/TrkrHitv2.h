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
  virtual ~TrkrHitv2() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrHitv2 class with adc = " << m_adc << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  // these set and get the energy before digitization
  virtual void addEnergy(const double edep);
  virtual double getEnergy();

  // after digitization, these are the adc values
  virtual void setAdc(const unsigned short adc) {m_adc = adc;}
  virtual unsigned int getAdc() ;

 protected:

  unsigned short m_adc = 0;
  ClassDef(TrkrHitv2, 1);
};

#endif //TRACKBASE_TRKRHITV2_H
