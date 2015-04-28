#ifndef __OBUFFER_H
#define __OBUFFER_H

#include "event_io.h"

#include "phenixTypes.h"
#include "oEvent.h"
#include "Event.h"

#ifdef WITHTHREADS
#include <pthread.h>
#endif


#ifndef __CINT__
class WINDOWSEXPORT oBuffer {
#else
class  oBuffer {
#endif

public:

  //** Constructors

  //  oBuffer( FILE *fpp, PHDWORD * , const int length
  //	 , const int iseq = 1, const int irun=1);

  oBuffer() {};

#ifndef WIN32
  oBuffer (int fd, PHDWORD * where, const int length,
		const int irun=1 , const int iseq=0 );
  
#ifdef WITHTHREADS

  oBuffer (int fd, const int length,
		const int irun=1 , const int iseq=0 );
#endif

#endif

  oBuffer (const char *filename, PHDWORD * where, const int length,
	   int &status,
	   const int irun=1 , const int iseq=0 );


  virtual  ~oBuffer();

  //  this creates a new event on the next address
  virtual int nextEvent( const int evtsize, const int etype =0, const int evtseq =0);


  // frame and packet adding

  virtual int addRawEvent(int *);

  virtual int addEvent(Event *);

  virtual int addFrame(PHDWORD *);

  virtual int  addPacket( const Packet *p);

  virtual int addUnstructPacketData(PHDWORD * data, 
		    const int length,
		    const int id,
		    const int wordsize,
		    const int hitformat);

  // now the transfer routine
  //int transfer ( dataProtocol * dprotocol );

  virtual int writeout ();

  // now the re-sizing of buffer
  virtual int setMaxSize( const int size);

  // and the query
  virtual int getMaxSize() const;

  // and the query`
  virtual unsigned long long getBytesWritten() const;

  virtual int addEoB();


protected:
  //  add end-of-buffer

  virtual int prepare_next();
#ifdef WITHTHREADS
  static void *writeThread(void * arg);
#endif
  typedef struct 
  { 
    int Length;
    int ID;
    int Bufseq;
    int Runnr;
    PHDWORD data[1];
  } *buffer_ptr;


  int we_are_threaded;
  buffer_ptr bptr;
  buffer_ptr bptr_being_written;

  buffer_ptr bptr0;
  buffer_ptr bptr1;
  PHDWORD *data_ptr;
  int current_index;
  int max_length;
  int max_size;
  unsigned int left;
  int sequence;
  int runnumber;
  unsigned long long byteswritten;
  oEvent *current_event;
  int eventsequence;
  int current_etype;
  int has_end;
  int dirty;
  int fd;
  int our_fd;
  int good_object;

#ifdef WITHTHREADS
#ifndef __CINT__
  pthread_t ThreadId;
#else
  int  ThreadId;
#endif
#endif
  int thread_arg[3];

};

#endif

