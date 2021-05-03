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
  TrkrHitv1(); 

  //! dtor
  virtual ~TrkrHitv1() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrHitV1 class with adc = " << m_adc  << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  virtual void addEnergy(const double edep) {m_edep += edep;}
  virtual double getEnergy() {return m_edep;}
  virtual void setAdc(const unsigned int adc) {m_adc = adc;}
  virtual unsigned int getAdc();

 protected:
  double m_edep = 0;
  unsigned int m_adc = 0;
  ClassDef(TrkrHitv1, 1);
};

#endif //TRACKBASE_TRKRHITV1_H
