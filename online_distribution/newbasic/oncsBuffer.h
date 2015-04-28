#ifndef __ONCSBUFFER_H
#define __ONCSBUFFER_H

#include "oncsEvent.h"
#include "BufferConstants.h"

#include "event_io.h"

#include <cstdio>

#ifndef __CINT__

class WINDOWSEXPORT oncsBuffer {
#else
class  oncsBuffer {
#endif

public:

  //** Constructors

  oncsBuffer( int *array, const int length);

  //  this creates a new event on the next address
  Event * getEvent();

  int buffer_swap();

  static int i4swap (const int in);
  static int i22swap (const int in);
  static short i2swap (const short in);

protected:
  typedef struct 
  { 
    int Length;
    int ID;
    int Bufseq;
    int Runnr;
    int data[1];
  } *buffer_ptr;

  buffer_ptr bptr;
  int *data_ptr;
  int buffer_size;
  int max_length;
  int current_index;
};

#endif
