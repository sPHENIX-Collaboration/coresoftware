// -*- c++ -*-
#ifndef __EVENTITERATOR_H__
#define __EVENTITERATOR_H__

#include "Event.h"


/**
   The Eventiterator creates and returns pointers to Event objects.
   This is the virtual base class of a variety of different iterators.
*/
#ifndef __CINT__
class WINDOWSEXPORT Eventiterator {
#else
class  Eventiterator {
#endif
public:
  /// The virtual destructor.
  inline virtual ~Eventiterator(){};

  virtual const char * getIdTag() const =0;

  virtual void identify(std::ostream& os = std::cout) const = 0;


/**
   this member function returns a pointer to the Event object, or
   NULL if there no events left.
*/   
  virtual Event *getNextEvent() = 0;


  /** the next member functions do nothing for most iterators
  except the etEventiterator. It's still easier to give the parent
  class the interfaces
  */

  virtual int *getNextEventData() {return 0;};
  virtual int releaseEventData() {return 0;};

  virtual void setBlockingMode(const int mode) {};
  virtual int getBlockingMode() const {return 0;};

  virtual void setSelectMode(const int mode) {};
  virtual int getSelectMode() const {return 0;};

  virtual void setSelectWords(const int i1, const int i2, const int i3, const int i4 ) {};
  virtual void getSelectWords(int val[]) const {};

  virtual const char * getCurrentFileName() const { return " "; };

  virtual int  setVerbosity(const int v) { return -1; }; // most iterators don't have the concept
  virtual int  getVerbosity() const { return -1; };


};

#endif /* __EVENTITERATOR_H__ */
