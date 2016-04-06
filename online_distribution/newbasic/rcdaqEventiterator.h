// -*- c++ -*-
#ifndef __RCDAQEVENTITERATOR_H__
#define __RCDAQEVENTITERATOR_H__

#include "Eventiterator.h"
#include "buffer.h"

#ifndef __CINT__
#include <string>

#include <arpa/inet.h>
#include <stdio.h>

#endif


#ifndef __CINT__
class WINDOWSEXPORT rcdaqEventiterator : public Eventiterator {
#else
class  rcdaqEventiterator : public Eventiterator {
#endif
public:

  virtual ~rcdaqEventiterator();
  rcdaqEventiterator();
  rcdaqEventiterator(const char *ip );
  rcdaqEventiterator(const char *ip, int &status);

  const char * getIdTag() const;
  virtual void identify(std::ostream& os = std::cout) const;


  Event *getNextEvent();

protected:
  int read_next_buffer();

  int setup(const char *ip, int &status);

  int readn (int fd, char *ptr, int nbytes);
  int writen (int fd, char *ptr, int nbytes);

  
  std::string _theIP;

  int _sockfd;
  PHDWORD initialbuffer[BUFFERSIZE];
  PHDWORD *bp;
  int allocatedsize;

  int current_index;
  int last_read_status;
  int buffer_size;
  buffer *bptr;
  struct sockaddr_in server;
  int _defunct;

};

#endif /* __RCDAQEVENTITERATOR_H__ */

