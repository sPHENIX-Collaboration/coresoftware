// -*- c++ -*-
#ifndef __TESTEVENTITERATOR_H__
#define __TESTEVENTITERATOR_H__

#include "Eventiterator.h"

#include "simpleRandom.h"

/** The testEventiterator creates Event objects with known properties
   from scratch. The events contain 1 frame and 3 packets. The first 2
   packets with id's 1001 and 1002 look like ADC's whose channel i
   contains the value i, which can be used to verify that a certain
   channel actually ends up where it should. The third packet 1003
   has four channels which give values according to gaussian distributions 
   with RMS 10,100,1000, and 10000, respectively. 
*/
#ifndef __CINT__
class WINDOWSEXPORT testEventiterator : public Eventiterator {
#else
class  testEventiterator : public Eventiterator {
#endif
public:
  /// The constructor.
  testEventiterator();
  /// The destructor.
  ~testEventiterator();

  const char * getIdTag() const;
  //#ifdef WIN32
  virtual void identify(std::ostream& os = std::cout) const;


/**
   this member function creates returns a pointer to the Event object
   (this iterator does not become empty).
*/   
  Event *getNextEvent();

private:
#ifndef WIN32
  simpleRandom *R;
#endif
  int current_event;

};

#endif /* __TESTEVENTITERATOR_H__ */

