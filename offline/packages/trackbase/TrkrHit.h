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

#include <climits>
#include <cmath>
#include <iostream>

/**
 * @brief Base class for hit object
 *
 * This is the empty virtual base class for a hit object.
 * Each subsystem should implement an inherited version
 * which contains the actual storage information.
 */
class TrkrHit : public PHObject
{
 public:

  //! dtor
  ~TrkrHit() override = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrHit base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  virtual void CopyFrom(const TrkrHit&)
  {}

  //! copy content from base class
  virtual void CopyFrom(TrkrHit*)
  {}

  // these set and get the energy before digitization
  virtual void addEnergy(const double) {}
  virtual double getEnergy() const { return 0; }

  // after digitization, these are the adc values
  virtual void setAdc(const unsigned int) {}
  virtual unsigned int getAdc() const { return 0; }
  /*
  virtual void setCrossing(const short int) {}
  virtual short int getCrossing() { return 0;}
  */

 protected:
  ClassDefOverride(TrkrHit, 1);
};

#endif  // TRACKBASE_TRKRHIT_H
