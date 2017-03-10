#ifndef __ONCSBUFFER_H
#define __ONCSBUFFER_H

#include "buffer.h"
#include "oncsEvent.h"
#include "BufferConstants.h"

#include <stdio.h>

#include <event_io.h>

#ifndef __CINT__

#define PRDFBUFFERID 0xffffffc0
#define ONCSBUFFERID 0xffffc0c0

class WINDOWSEXPORT oncsBuffer : public buffer {
#else
  class  oncsBuffer : public buffer {
#endif

public:

  //** Constructors

  oncsBuffer( PHDWORD *array, const PHDWORD length);

  //  this creates a new event on the next address
  Event * getEvent();

  int buffer_swap();

  int *getEventData() { return 0;};
  int isGood() const { return 1;};

  static int i4swap (const int in);
  static int i22swap (const int in);
  static short i2swap (const short in);

protected:
  typedef struct 
  { 
    int Length;
    unsigned int ID;
    int Bufseq;
    int Runnr;
    int data[];
  } *buffer_ptr;

  buffer_ptr bptr;
  int *data_ptr;
  int buffer_size;
  int max_length;
  int current_index;
};

#endif
