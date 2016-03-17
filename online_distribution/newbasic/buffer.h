#ifndef __BUFFER_H
#define __BUFFER_H

#include "phenixTypes.h"
#include "BufferConstants.h"
#include "Event.h"



#ifndef __CINT__
class WINDOWSEXPORT buffer {
#else
class  buffer {
#endif

public:

  //** Constructors

  buffer();
  virtual ~buffer();

  //  this creates a new event on the next address
  virtual Event * getEvent() = 0;

  virtual int * getEventData() = 0;

  virtual int isGood() const = 0;

  static int makeBuffer( PHDWORD *bp, const int allocatedsize, buffer **bptr);
  static int i4swap (const int in);
  static unsigned int u4swap (const unsigned int in);
  static int i22swap (const int in);
  static short i2swap (const short in);

protected:
};

#endif
