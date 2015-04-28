#ifndef __OAMLBUFFER_H__
#define __OAMLBUFFER_H__

#include "olzoBuffer.h"

#define CTRL_BEGINRUN 1
#define CTRL_ENDRUN   2
#define CTRL_DATA     3
#define CTRL_CLOSE    4


#ifndef __CINT__
class WINDOWSEXPORT oamlBuffer : public olzoBuffer{
#else
class  oamlBuffer : public oBuffer{
#endif

public:

  //** Constructors

  oamlBuffer (const char *host, const int port, PHDWORD * where, 
	     const int length,
	     const int irun=0, 
	     const int iseq=0 );

  // this one takes the form hostname:port
  
  oamlBuffer (const char *hostport,  PHDWORD * where, 
	     const int length,
	     const int irun=0, 
	     const int iseq=0 );


  virtual  ~oamlBuffer();


  virtual int writeout ();


protected:

  int connect_aml();
  int begin_run();

  int has_begun;
  int RunNumber;
  int DontOverrideRunNumber;

  static int readn(int , char *, int);
  static int writen(int , char *, int);


  char HostName[256];
  int ThePort;

  int sockfd;


};



#endif /* __OAMLBUFFER_H__ */

