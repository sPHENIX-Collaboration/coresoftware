
#include <remote_msg_buffer.h>

// in the constructor we allocate the specified amount of bytes for
// the output array.

remote_msg_buffer::remote_msg_buffer(const char *host, const int port, const int msglen)
  : msg_buffer(msglen)
{
  // in the constructor we replace COUT's streambuf with this one.
  // we remember the old one because we may want to replace it.

#ifdef RDBUF_ACCEPTS_STREAMBUF
  original_streambuf = COUT.rdbuf ( (STREAMBUF *) this);
#else
  original_streambuf = COUT.rdbuf ();
  COUT = (STREAMBUF *) this;
#endif

  ThePort = port;
  TheHost = new char [ strlen(host) +1 ];
  strcpy (TheHost, host);

  p_host = gethostbyname(TheHost);
  COUT << p_host->h_name << ENDL;

  memset ( (void*) &server_addr, 0, sizeof(server_addr) );
  server_addr.sin_family = AF_INET;

  //bcopy(p_host->h_addr, &(server_addr.sin_addr.s_addr), p_host->h_length);
  memcpy (&(server_addr.sin_addr.s_addr), p_host->h_addr, p_host->h_length);

  server_addr.sin_port = htons(ThePort);

}

// in the destructor, we restore the original streambuf for COUT.
remote_msg_buffer::~remote_msg_buffer()
{
#ifdef RDBUF_ACCEPTS_STREAMBUF
  COUT.rdbuf ( original_streambuf);
#else
  COUT= original_streambuf;
#endif
  delete [] TheHost;

}

// this is the function which is called as a result of the << std::endl 
// and flush() operations
#if defined(LVL2_WINNT) || defined(STREAMBUF_NEW_IOSTREAM)
int remote_msg_buffer::pubsync ()
#else
int remote_msg_buffer::sync ()
#endif
{ 

  msgProfile mp;
  int length, ix;


  char *x = format(&length, &mp);
  int ret = 0;

#if defined(LVL2_WINNT) || defined(STREAMBUF_NEW_IOSTREAM)
  for (ix =0; ix < length; ix++)
    {
      original_streambuf->sputc( x[ix] );
    }
  ret = original_streambuf->pubsync();
#else
  for (ix =0; ix < length; ix++)
    {
      original_streambuf->overflow( x[ix] );
    }
  ret = original_streambuf->sync();
#endif

  send_string(x, length);

  delete [] x;

  pos = 0;
  
  return ret;
}

void remote_msg_buffer::send_string(const char *x, const int len)
{
  

  if ( (sockfd = socket(AF_INET, SOCK_STREAM, 0) ) < 0 )
    {
      return;
    }

  if ( connect(sockfd, (struct sockaddr*) &server_addr
	       , sizeof(server_addr)) < 0 ) 
    {
      return;
    }
  writen (sockfd, x, len);

  close (sockfd);
}

int remote_msg_buffer::writen (int fd, const char *ptr, int nbytes)
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
