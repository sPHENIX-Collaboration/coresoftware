#ifndef __PRDFBUFFER_H
#define __PRDFBUFFER_H

#include <buffer.h>
#include <phenixTypes.h>
#include <BufferConstants.h>
#include <Event.h>



#ifndef __CINT__
class WINDOWSEXPORT prdfBuffer : public buffer {
#else
  class  prdfBuffer : public buffer  {
#endif

public:

  //** Constructors

  prdfBuffer();
  prdfBuffer( PHDWORD *array, const int length);
   ~prdfBuffer();

  //  this creates a new event on the next address
  Event * getEvent();

  int * getEventData();

  int isGood() const { return is_good; } ;

  int buffer_swap();
  int frame_swap(PHDWORD * fp, const int eventlength);

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
