#include "A_Event.h"
#include "packet.h"
#include "Cframe.h"
#include "framePackets.h"
#include "dataBlock.h"

#include "frameRoutines.h"
#include "frameHdr.h"

#include "fileEventiterator.h"
#include "testEventiterator.h"
#include "phenixTypes.h"
#include "oBuffer.h"
#include "oEvent.h"


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>



#include <stdlib.h>
#include <signal.h>
#include <unistd.h>

#include <stdio.h>

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define DDEVENTITERATOR 1
#define FILEEVENTITERATOR 2
#define TESTEVENTITERATOR 3
#define DDPOOL 4
#define DFILE 5
#define DNULL 6




class X_Event : public A_Event 
{

public:
  // constructors and destructors
  X_Event( PHDWORD *);
  X_Event( int *);
  ~X_Event();

   // info & debug utils
   
   int change_id(const int oldid, const int newid);


};

X_Event::X_Event (PHDWORD *data)
  : A_Event(data)
{ }

X_Event::X_Event (int *data)
  : A_Event(data)
{ }

X_Event::~X_Event ()
{ }

int X_Event::change_id (const int oldid, const int newid)
{
  
  int i = 0;
  PHDWORD *fp;
  PHDWORD *pp;
  UINT ids = oldid;
  
  while ( fp = framelist[i++] )
    {
      if ( ( pp = findFramePacketId (fp, ids) ) !=  ptrFailure) 
	{
	  setPacketId (pp, newid);
	}
    }
  return 0;

}




#if defined(SunOS) || defined(Linux) || defined(OSF1)
void sig_handler(int);
#else
void sig_handler(...);
#endif


void exitmsg()
{
  COUT << "** usage: changeid infile outfile id1 id2 " << std::endl;
  exit(0);
}

void exithelp()
{
  COUT <<  std::endl;

}

// The global pointer to the Eventiterator (we must be able to 
// get at it in the signal handler) 
Eventiterator *it;



int 
main(int argc, char *argv[])
{
  int c;
  int status = 0;

  int eventnr = 0;

  extern char *optarg;
  extern int optind;

  PHDWORD  *buffer;
  oBuffer *ob;
  int fd;
  int buffer_size = 2000000;

  if (argc < 2) exitmsg();

  while ((c = getopt(argc, argv, "s:d:n:w:vhi")) != EOF)
    {
      switch (c) 
	{
	  // the -s (source type) switch
	case 'h':
	  exithelp();
	  break;

	default:
	  break;
	}
    }



  // install some handlers for the most common signals
  signal(SIGKILL, sig_handler);
  signal(SIGTERM, sig_handler);
  signal(SIGINT,  sig_handler);

  // see if we can open the file


  it = new fileEventiterator(argv[optind], status);
       

  if (status)
    {
      delete it;
      COUT << "Could not open input stream" << std::endl;
      exit(1);
    }

  unlink ( argv[optind+1] );
  fd = open(argv[optind+1], O_WRONLY | O_CREAT | O_EXCL | O_LARGEFILE , 
		  S_IRWXU | S_IROTH | S_IRGRP );
  if ( fd < 0) 
    {
      COUT << "Could not open file: " <<  argv[optind+1] << std::endl;
      exit (1);
    }

  buffer = new PHDWORD [buffer_size];
  
  ob = new oBuffer (fd, buffer, buffer_size);


  int paircount = 0;
  int idold[1000];
  int idnew[1000];

  int argind = 0;
  while ( optind +2 + argind +1 < argc) 
    {
      !sscanf(argv[optind +2 + argind]    , "%d", &idold[paircount]); 
      !sscanf(argv[optind +2 + argind + 1], "%d", &idnew[paircount]); 
      COUT << "changing  " << idold[paircount] << " -> " << idnew[paircount] << std::endl;
      argind+=2;
      paircount++;
    }

	

  // ok. now go through the events

  Event *evt;
  X_Event *nevt;

  int *eb;
  int nw;
  int j;
  while (  ( evt = it->getNextEvent())  )
    {


      eb = new int [evt->getEvtLength() +100];
      evt->Copy (eb, evt->getEvtLength() +100, &nw);
      nevt = new X_Event(eb);
      for (j=0; j<paircount; j++)
	nevt->change_id(idold[j] , idnew[j]);

      ob->addEvent(nevt);


      eventnr++;

      delete evt;
      delete nevt;
      delete [] eb;

    }
  delete it;


  delete ob;
  close(fd);
      
  return 0;
}
	

#if defined(SunOS) || defined(Linux) || defined(OSF1)
void sig_handler(int i)
#else
  void sig_handler(...)
#endif
{
  COUT << "sig_handler: signal seen " << std::endl;
  if (it) delete it;
  exit(0);
}


  
