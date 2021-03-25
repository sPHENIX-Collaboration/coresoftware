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
  TrkrHit(); 

  //! dtor
  virtual ~TrkrHit() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrHit base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  // these set and get the energy before digitization
  void addEnergy(const double edep) {
    m_adc += (unsigned int) (edep* TrkrDefs::EdepScaleFactor); 
  }
  double getEnergy() {
    return ((double) m_adc)  / TrkrDefs::EdepScaleFactor;
  }
  // after digitization, these are the adc values
  void setAdc(const unsigned int adc) {m_adc = adc;}
  unsigned int getAdc() { return m_adc;}

 protected:
  //double m_edep;
  unsigned int m_adc;
  ClassDef(TrkrHit, 1);
};

#endif //TRACKBASE_TRKRHIT_H
