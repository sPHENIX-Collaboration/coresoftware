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
  buffer( PHDWORD *array, const int length);
  virtual ~buffer();

  //  this creates a new event on the next address
  virtual Event * getEvent();

  virtual int * getEventData();

  virtual int isGood() const { return is_good; } ;

  int buffer_swap();
  int frame_swap(PHDWORD * fp, const int eventlength);

  static int i4swap (const int in);
  static unsigned int u4swap (const unsigned int in);
  static int i22swap (const int in);
  static short i2swap (const short in);

protected:
  typedef struct 
  { 
    unsigned int Length;
    unsigned int ID;
    int Bufseq;
    int Runnr;
    PHDWORD data[1];
  } *buffer_ptr;

  buffer_ptr bptr;
  PHDWORD *data_ptr;
  int buffer_size;
  int max_length;
  int current_index;
  int is_good;
};

#endif
