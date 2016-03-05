
#ifndef __MSG_BUFFER_H__
#define __MSG_BUFFER_H__


#include "event_io.h"

#ifdef __CINT__
#include <iostream>
#endif


#ifndef __CINT__

#include <cstdio>
#include <cstdarg>
#include <string>
#include <ctime>
#include <sys/timeb.h>

#endif


#include "msg_control.h"

// if we find that our buffer is too small, we extend it
// by this amount.
#define MSG_EXTENSION_AMOUNT 32


typedef struct msgProfile {
  int type;
  int source;
  int severity;
  int reserved1;
  int reserved2;
  int reserved3;
  char sourcecomponent[256];
} msgProfile;



//class msg_control;


/** 
This is the base class for the type of message buffers which we want to 
use for messaging. It can easily be subclassed. 


This base class defines the "static" behavior of all the msg\_buffer
classes derived from it. It provides the parsing routine which decodes
the profile tags, provides a "add\_date" routine to conveniently
prepend a date tag to each message if so desired, and deals with the
msg\_control objects in the system which it is a friend of. If created, 
the msg\_buffer object tells all msg\_control object to go active and produce
the profiling information; otherwise they maintain silence. 

This class is meant to be subclassed, and the dispatch mechanism which
you want has to be provided by the subclass.

*/

#ifdef __CINT__
class msg_buffer : public streambuf {
#else
class msg_buffer : public STREAMBUF {
#endif

protected:
  char *oBuffer;
  int pos;
  int maximum_position;
  msg_control *m;

  virtual char * format(int * length, msgProfile *mp);

  

  public:
    msg_buffer (const int msglen=256);
    virtual ~msg_buffer();

    virtual int overflow (int ch);
    virtual int sync ();
};

#ifndef __CINT__
void streambuf_add_date (STREAMBUF * sb);
#endif


#endif /* __MSG_BUFFER_H__ */
