// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_EVENTHEADERV2_H
#define FFAOBJECTS_EVENTHEADERV2_H

/*!
 * \file EventHeaderv2.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "EventHeaderv1.h"

#include <ctime>         // for time_t
#include <iostream>       // for cout, ostream

class PHObject;

//! simple event header with ID and time
class EventHeaderv2: public EventHeaderv1
{
 public:

  //! ctor
  EventHeaderv2();

  //! dtor
  virtual ~EventHeaderv2() = default;

  //! clone
  PHObject *CloneMe() const 
  { return new EventHeaderv2(*this); }

  ///  Clear Event
  void Reset();

  /*!
   * identify Function from PHObject
   * @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const;

  //! bunch crossing
  void set_BunchCrossing( int64_t value ) 
  { m_bunchCrossing = value; }
  
  //! bunch crossing
  int64_t get_BunchCrossing() const
  { return m_bunchCrossing; }
  
  private: 

  //! bunch crossing id
  int64_t m_bunchCrossing = 0;

  ClassDef(EventHeaderv2,1)
};

#endif
