// -*- c++ -*-
#ifndef __RCDAQEVENTITERATOR_H__
#define __RCDAQEVENTITERATOR_H__

#include "Eventiterator.h"
#include "oncsBuffer.h"

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
  rcdaqEventiterator(const char *ip = "127.0.0.1");
  rcdaqEventiterator(const char *ip, int &status);

  char * getIdTag() const;
  virtual void identify(std::ostream& os = std::cout) const;


  Event *getNextEvent();

protected:
  int read_next_buffer();

  int setup(const char *ip, int &status);
  
  std::string _theIP;

  int _sockfd;
  int initialbuffer[BUFFERSIZE];
  int *bp;
  int allocatedsize;

  int current_index;
  int last_read_status;
  int buffer_size;
  oncsBuffer *bptr;
  struct sockaddr_in server;

};

#endif /* __RCDAQEVENTITERATOR_H__ */

