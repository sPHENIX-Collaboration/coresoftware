// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_EVENTHEADERV2_H
#define FFAOBJECTS_EVENTHEADERV2_H

/*!
 * \file EventHeaderv2.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "EventHeaderv1.h"

#include <cstdint>   // for int64_t
#include <iostream>  // for cout, ostream

class PHObject;

//! simple event header with ID and time
class EventHeaderv2 : public EventHeaderv1
{
 public:
  //! ctor
  EventHeaderv2() = default;

  //! dtor
  virtual ~EventHeaderv2() = default;

  //! clone
  PHObject* CloneMe() const override
  {
    return new EventHeaderv2(*this);
  }

  ///  Clear Event
  void Reset() override;

  /*!
   * identify Function from PHObject
   * @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const override;

  //! bunch crossing
  void set_BunchCrossing(int64_t value) override
  {
    m_bunchCrossing = value;
  }

  //! bunch crossing
  int64_t get_BunchCrossing() const override
  {
    return m_bunchCrossing;
  }

 private:
  //! bunch crossing id
  int64_t m_bunchCrossing = 0;
// rootcling and clang complain about inconsistent overrides in the ClassDef
// this can be supressed with ignoring -Winconsistent-missing-override
// this pragma is not known to gcc, so we need an #ifdef __clang__ here
#pragma GCC diagnostic push
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#endif
  ClassDef(EventHeaderv2, 2)
#pragma GCC diagnostic pop
};

#endif
