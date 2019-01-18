/**
 * @file trackbase/TrkrHit.h
 * @author D. McGlinchey
 * @date 4 June 2018
 * @brief Base class for hit object
 */
#ifndef TRACKBASE_TRKRHIT_H
#define TRACKBASE_TRKRHIT_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>

/**
 * @brief Base class for hit object
 *
 * This is the empyt virtual base class for a hit object.
 * Each subsystem should implement an inherited version
 * which contains the actual storage information.
 */
class TrkrHit : public PHObject
{
 public:
  //! ctor
  TrkrHit() {}
  //! dtor
  virtual ~TrkrHit() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrHit base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  void addEnergy(const double edep) {Edep += edep;}
  double getEnergy() {return Edep;}
  void setAdc(const unsigned int adc) {Adc = adc;}
  unsigned int getAdc() {return Adc;}

 protected:
  double Edep;
  unsigned int Adc;

  ClassDef(TrkrHit, 1);
};

#endif //TRACKBASE_TRKRHIT_H
