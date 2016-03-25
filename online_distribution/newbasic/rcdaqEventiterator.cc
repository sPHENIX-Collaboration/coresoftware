//
// rcdaqeventIterator   mlp 4/19/1997
//
// this iterator reads events froma data file. 

#define MONITORINGPORT 9930

#include "rcdaqEventiterator.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#include "oncsBuffer.h"
#include "gzbuffer.h"
#include "lzobuffer.h"
#include "Event.h"

#include <stddef.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


#include <netdb.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/socket.h>


using namespace std;


// there are two similar constructors, one with just the
// filename, the other with an additional status value
// which is non-zero on return if anything goes wrong. 

rcdaqEventiterator::~rcdaqEventiterator()
{
     if (_sockfd) close (_sockfd);
     if (bp != NULL ) delete [] bp;
     if (bptr != NULL ) delete bptr;
}  


rcdaqEventiterator::rcdaqEventiterator()
{
  string host = "localhost";
    
  if ( getenv("RCDAQHOST")  )
    {
      host = getenv("RCDAQHOST");
    }

  int status;
  setup (host.c_str(), status);
}  

rcdaqEventiterator::rcdaqEventiterator(const char *ip)
{
  int status;
  setup (ip, status);
}  

rcdaqEventiterator::rcdaqEventiterator(const char *ip, int &status)
{
  setup (ip, status);
}  


int rcdaqEventiterator::setup(const char *ip, int &status)
{
  _defunct = 0;

  
  struct hostent *p_host;
  p_host = gethostbyname(ip);
  
  //  std::cout << __FILE__ << " " << __LINE__ << "  " << ip << "  " << p_host->h_name << "  " << p_host->h_addr << std::endl;

  
  if ( ! p_host ) 
    {
      status = -2;
      _defunct = 1;
      return -2;
    }

  //  std::cout << p_host->h_name << "  " << p_host->h_addr << std::endl;

  bptr = 0;
  bp = 0;
  allocatedsize = 0;
  _theIP = p_host->h_name;
  status = 0;
  last_read_status = 0;
  current_index = 0;

  memset((char *) &server, 0, sizeof(server));
  server.sin_family = AF_INET;
  bcopy(p_host->h_addr, &(server.sin_addr.s_addr), p_host->h_length);
  server.sin_port = htons(MONITORINGPORT);

  return 0;
}  

void rcdaqEventiterator::identify (OSTREAM &os) const
{ 
  os << getIdTag();
  if ( _defunct ) os << " *** defunct";
  os << std::endl;

};

const char * rcdaqEventiterator::getIdTag () const
{ 
  static char line[180];
  strcpy (line, " -- rcdaqEventiterator reading from ");
  strcat (line, _theIP.c_str());
  return line;
};


// and, finally, the only non-constructor member function to
// retrieve events from the iterator.

Event * rcdaqEventiterator::getNextEvent()
{
  if ( _defunct ) return 0;
  Event *evt = 0;

  // if we had a read error before, we just return
  if (last_read_status) return NULL;

  // see if we have a buffer to read
  if (bptr == 0) 
    {
      if ( (last_read_status = read_next_buffer()) !=0 )
	{
	  return NULL;
	}
    }

  while (last_read_status == 0)
    {
      if (bptr) evt =  bptr->getEvent();
      if (evt) return evt;

      last_read_status = read_next_buffer();
    }

  return NULL;

}

// -----------------------------------------------------
// this is a private function to read the next buffer
// if needed. 

int rcdaqEventiterator::read_next_buffer()
{
  int ip = 0;
  if (bptr) 
    {
      delete bptr;
      bptr = 0;
    }

  _sockfd  = socket(AF_INET, SOCK_STREAM, 0);
  if ( _sockfd < 0) return 0;

  if ( connect(_sockfd, (struct sockaddr*) &server, sizeof(server)) < 0 ) 
	{
	  //std::cout << "error in connect" << std::endl;
	  close (_sockfd);
	  usleep(1000);  // we just slow down a bit to limit the rate or retries
	  return 0;
	}
  
  // say that this is our max size
  int flag = htonl(64*1024*1024);
  
  int status = writen (_sockfd,(char *)  &flag, 4);
  if ( status < 0)
    {
      close (_sockfd);
      return 0;
    }

  
  int sizetobesent;
  status = readn (_sockfd, (char *) &sizetobesent, 4);
  if ( status < 0)
    {
      close (_sockfd);
      return 0;
    }

  buffer_size = ntohl(sizetobesent);
  int i;
  if (bp) 
    {
      if  (buffer_size > allocatedsize*4)
	{
	  delete [] bp;
	  i = (buffer_size +8191) /8192;
	  allocatedsize = i * 2048;
	  bp = new PHDWORD[allocatedsize];
	}
    }
  else
    {
      i = (buffer_size +8191) /8192;
      allocatedsize = i * BUFFERBLOCKSIZE/4;
      bp = new PHDWORD[allocatedsize];

    }

  status = readn ( _sockfd, (char *) bp, buffer_size);
  if ( status < 0)
    {
      close (_sockfd);
      return 0;
    }

  int ackvalue = htonl(101);
  writen (_sockfd,(char *)  &ackvalue, 4);
  close (_sockfd);
  
  return buffer::makeBuffer( bp, allocatedsize, &bptr);

}

int rcdaqEventiterator::readn (int fd, char *ptr, int nbytes)
{
  int nleft, nread;
  nleft = nbytes;
  while ( nleft>0 )
    {
      nread = recv (fd, ptr, nleft, MSG_NOSIGNAL);
      if ( nread < 0 )
	{
	  return nread;
	}
      else if (nread == 0) 
	break;
      nleft -= nread;
      ptr += nread;
    }
  return (nbytes-nleft);
}


int rcdaqEventiterator::writen (int fd, char *ptr, int nbytes)
{
  int nleft, nwritten;
  nleft = nbytes;
  while ( nleft>0 )
    {
      nwritten = send (fd, ptr, nleft, MSG_NOSIGNAL);
      if ( nwritten < 0 )
	{
	  return nwritten;
	}
      nleft -= nwritten;
      ptr += nwritten;
    }
  return (nbytes-nleft);
}
