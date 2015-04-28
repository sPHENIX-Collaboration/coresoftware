
#include "oamlBuffer.h"
#include "BufferConstants.h"

#include <unistd.h>
#include <netdb.h>
#include <cstring>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>

#if defined(SunOS) || defined(OSF1)
#include <strings.h>
#endif


// the constructor first ----------------
oamlBuffer::oamlBuffer ( const char * host, const int port, PHDWORD * where, 
		      const int length, 
		      const int irun, 
		      const int iseq): 
  olzoBuffer(0,where,length,irun,iseq)
{
  good_object = 1;

  has_begun = 0;  // we will defer the connection until we have the first buffer to write


  ThePort = port;
  strcpy(HostName, host);

  if ( connect_aml() ) good_object =0;

  DontOverrideRunNumber = 0;
  RunNumber = irun;
  if ( irun) 
    {
      DontOverrideRunNumber = 1;
    }
}


oamlBuffer::oamlBuffer ( const char * hostport,  PHDWORD * where, 
		      const int length, 
		      const int irun, 
		      const int iseq): 
  olzoBuffer(0,where,length,irun,iseq)
{
  good_object = 1;

  has_begun = 0;  // we will defer the connection until we have the first buffer to write

  char *t = new char[strlen(hostport)+1];
  strcpy ( t, hostport);

  char *token = strtok(t, ":");
  if (token)
    {
      strcpy( HostName, token);
      token = strtok(0, ":");
      if (token)
	{
	  sscanf(token, "%d", &ThePort);
	}
    }

  delete [] t;

  if ( connect_aml() ) good_object =0;

  DontOverrideRunNumber = 0;
  RunNumber = irun;
  if ( irun) 
    {
      DontOverrideRunNumber = 1;
    }

}


int oamlBuffer::connect_aml ()
{

  struct sockaddr_in server_addr;
  struct hostent *p_host;
  p_host = gethostbyname(HostName);

  if ( ! p_host) return -3;

  bzero( (char*) &server_addr, sizeof(server_addr) );
  server_addr.sin_family = AF_INET;
  bcopy(p_host->h_addr, &(server_addr.sin_addr.s_addr), p_host->h_length);
  server_addr.sin_port = htons(ThePort);


  if ( (sockfd = socket(AF_INET, SOCK_STREAM, 0) ) < 0 )
    {
      std::cout << "error in socket" << std::endl;
      good_object = 0;
      return -1;
    }

  int xs = 512*1024;
  
  int s = setsockopt(sockfd, SOL_SOCKET, SO_SNDBUF,
		     &xs, 4);
  if (s) std::cout << "setsockopt status = " << s << std::endl;

  if ( connect(sockfd, (struct sockaddr*) &server_addr
	       , sizeof(server_addr)) < 0 ) 
    {
      std::cout << "error in connect" << std::endl;
      good_object = 0;
      return -1;
    }

  return 0;
}

int oamlBuffer::begin_run ()
{
  int i;
  int controlword = htonl(CTRL_BEGINRUN);
  int hrun;
  if ( DontOverrideRunNumber )
    {
      hrun = htonl(RunNumber);  // this is our own, cached value from the constructor
    }
  else
    {
      hrun = htonl(runnumber);  // this is the one which the oBuffer parent maintains 
    }
  writen (sockfd,(char *)  &controlword, 4);
  writen (sockfd,(char *) &hrun, 4);
  readn (sockfd, (char *) &i, 4);
  has_begun = 1;

  return 0;

}

// ----------------------------------------------------------
// returns the number of bytes written, including record wasted space.
//
int oamlBuffer::writeout()
{

  int cstatus = 0;

  if (! dirty) return 0;

  if ( ! good_object) return 1;
  if ( ! has_begun ) cstatus = begin_run();

  if ( cstatus ) return cstatus;

  if (! has_end) addEoB();

  lzo_uint outputlength_in_bytes = outputarraylength*4-16;
  lzo_uint in_len = bptr->Length; 

  lzo1x_1_12_compress( (lzo_byte *) bptr,
		       in_len,  
		       (lzo_byte *)&outputarray[4],
		       &outputlength_in_bytes,wrkmem);


  outputarray[0] = outputlength_in_bytes +4*BUFFERHEADERLENGTH;
  outputarray[1] =  LZO1XBUFFERMARKER;
  outputarray[2] = bptr->Bufseq;
  outputarray[3] = bptr->Length;


  int controlword = htonl(CTRL_DATA);
  //  std::cout << __LINE__ << "  sent control word" << std::endl;
  writen (sockfd,(char *)  &controlword, 4);
  int i = htonl(outputarray[0]);
  writen (sockfd,(char *) &i, 4);
  writen (sockfd,(char *)outputarray, outputarray[0] );
      
  //  std::cout << __LINE__ << "  waiting for ack" << std::endl;
  readn (sockfd, (char *) &i, 4);
  //  std::cout << __LINE__ << "  got it" << std::endl;
    


  dirty = 0;
  byteswritten += bptr->Length;
  return 0;
}


// ----------------------------------------------------------
oamlBuffer::~oamlBuffer()
{
  int wstatus =   writeout();

  if ( ! wstatus) 
    {
      
      int i;
      int controlword = htonl(CTRL_ENDRUN);
      writen (sockfd, (char *)&controlword, 4);
      
      
      readn (sockfd, (char *) &i, 4);
      
      controlword = htonl(CTRL_CLOSE);
      writen (sockfd, (char *)&controlword, 4);
    }
  dirty = 0;
  close (sockfd);
      
}

int oamlBuffer::readn (int fd, char *ptr, int nbytes)
{
  int nleft, nread;
  nleft = nbytes;
  while ( nleft>0 )
    {
      nread = read (fd, ptr, nleft);
      if ( nread < 0 ) 
	return nread;
      else if (nread == 0) 
	break;
      nleft -= nread;
      ptr += nread;
    }
  return (nbytes-nleft);
}


int oamlBuffer::writen (int fd, char *ptr, int nbytes)
{
  int nleft, nwritten;
  nleft = nbytes;
  while ( nleft>0 )
    {
      nwritten = write (fd, ptr, nleft);
      if ( nwritten < 0 ) 
	return nwritten;

      nleft -= nwritten;
      ptr += nwritten;
    }
  return (nbytes-nleft);
}
