//
// rcdaqeventIterator   mlp 4/19/1997
//
// this iterator reads events froma data file. 

#define MONITORINGPORT 9930

#include "rcdaqEventiterator.h"
#include <stdio.h>
#include <iostream>
#include "oncsEvent.h"
#include <stddef.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

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
  _sockfd  = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
  bptr = 0;
  bp = 0;
  allocatedsize = 0;
  if (_sockfd >0 ) 
    {
      _theIP = ip;
      status = 0;
      last_read_status = 0;
      current_index = 0;
    }
  else
    {
      status = 1;
      last_read_status = 1;
    }

  memset((char *) &server, 0, sizeof(server));
  server.sin_family = AF_INET;
  if (inet_aton(_theIP.c_str(), &server.sin_addr)==0) 
    {
      fprintf(stderr, "inet_aton() failed\n");
    }

  server.sin_port = htons(MONITORINGPORT);
  return 0;
}  

void  
rcdaqEventiterator::identify (OSTREAM &os) const
{ 
  os << getIdTag() << std::endl;

};

char * rcdaqEventiterator::getIdTag () const
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


  // say that this is our max size
  int flag = 64*1024*1024;
  socklen_t slen=sizeof(server);
  if (sendto(_sockfd, &flag, sizeof(int), 0, (sockaddr *) &server, slen)==-1)
    {
      perror("sendto 1 " );
    }
      


  int xc;
  int total;

  //receive the total number of bytes we are going to get 
  slen=sizeof(server);
  if ( ( xc = recvfrom(_sockfd, &total , 4 , 0, (sockaddr *) &server, &slen) ) ==-1)
    {
      perror("receive total number words " );
    }
  buffer_size = ntohl(total);


  // read the first record into initialbuffer
      
  //  COUT << "reading next buffer buffersize = " << buffer_size << std::endl; 


  int i;
  if (bp) 
    {
      if  (buffer_size > allocatedsize*4)
	{
	  delete [] bp;
	  i = (buffer_size +8191) /8192;
	  allocatedsize = i * 2048;
	  bp = new int[allocatedsize];
	}
    }
  else
    {
      i = (buffer_size +8191) /8192;
      allocatedsize = i * 2048;
      bp = new int[allocatedsize];
    }

  int max_sock_length = 48*1024;

  char *cp = (char *) bp;

  // now we read records until the whole buffer is read 
  while ( ip < buffer_size)
    {
      // read the next record
      
      
      slen=sizeof(server);
      if ( ( xc = recvfrom(_sockfd, cp , max_sock_length , 0, (sockaddr *) &server, &slen) ) ==-1)
	{
	  perror("receive " );
	}
      cp+=xc;
      ip+=xc; 
      //      cout << "sending ack " << ip  << endl;
      int ackvalue = 101;
      if (sendto(_sockfd, &ackvalue , sizeof(int), 0, (sockaddr *) &server, slen)==-1)
      	{
      	  perror("sendto 1 " );
      	}
      

    }

  //  cout << "received: " << ip << "  bytes" << endl;
  // and initialize the current_index to be the first event
  bptr = new oncsBuffer ( bp, allocatedsize );
  return 0;
}

