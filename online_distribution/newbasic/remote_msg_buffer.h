#ifndef __REMOTE_MSG_BUFFER_H__
#define __REMOTE_MSG_BUFFER_H__



#include "msg_buffer.h"

#ifndef __CINT__

#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/times.h>

#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#endif



/** This is the "filter" msg\_buffer class which allows you to filter
messages based on their profile. Its default behavior is to let all messages pass. 
You can use the member functions to tailor the kind of messages filtered and passed on.

*/

class remote_msg_buffer : public   msg_buffer {

protected:

#ifndef __CINT__
  STREAMBUF * original_streambuf;
  int ThePort;
  char *TheHost;

  int sockfd;
  struct sockaddr_in server_addr;
  struct hostent *p_host;

  void send_string(const char *x, const int len);
  int writen (int fd, const char *ptr, int nbytes);
#endif


public:


  /** The msglen parameter specifies the initial length of the 
      message string which is kept internally. If you exceed the length, 
      it is automautically extended.
  */
  remote_msg_buffer (const char *host = "va010.phenix.bnl.gov",
		     const int port = 8400, const int msglen=256);

  /// the virtual destructor
  virtual ~remote_msg_buffer();

  /// the sync function overrides the streambuf's sync function

  // mlp -- the pubsync is what's needed for the new
  // iostream libraries - sync will no longer be a public method. 
  // for now we leave it as it was.

#ifdef STREAMBUF_NEW_IOSTREAM
  virtual int pubsync ();
#else
  virtual int sync ();
#endif


};


#endif /* __REMOTE_MSG_BUFFER_H__ */
