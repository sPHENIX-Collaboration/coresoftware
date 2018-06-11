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
 * This is the base class for a hit object.
 * Each subsystem should implement an inherited version
 * which contains the actual storage information.
 */
class TrkrHit : public PHObject
{
 public:
  //! ctor
  TrkrHit() { m_key = TrkrDefs::HITKEYMAX; }
  //! dtor
  virtual ~TrkrHit() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrHit base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  /**
   * @brief Set the key for this hit
   * @param[in] key
   */
  void setKey(const TrkrDefs::hitkey key) { m_key = key; }

  /**
   * @brief Get the key for this hit
   * @param[out] key
   */
  TrkrDefs::hitkey getKey() const { return m_key; }

 protected:
  TrkrDefs::hitkey m_key;
  ClassDef(TrkrHit, 1);
};

#endif //TRACKBASE_TRKRHIT_H
