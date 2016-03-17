#include "fileEventiterator.h"
#include "testEventiterator.h"
#include "phenixTypes.h"
#include "A_Event.h"
#include "ogzBuffer.h"
#include "EventTypes.h"
#include "oEvent.h"

#include <cstdlib>
#include <unistd.h>
#include <cstdio>

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define DDEVENTITERATOR 1
#define FILEEVENTITERATOR 2
#define TESTEVENTITERATOR 3
#define DDPOOL 4
#define DFILE 5
#define DNULL 6

#if defined(SunOS) || defined(Linux) || defined(OSF1)
#include <cstring>
#else
#include <bstring.h>
#endif



#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#if defined(SunOS) || defined(Linux) || defined(OSF1)
void sig_handler(int);
#else
void sig_handler(...);
#endif


void exitmsg()
{
  COUT << "** usage: eventcombiner outputfile inputfile1 inputfile2 ..." << std::endl;
  COUT << "          eventcombiner -h for more help" << std::endl;
  exit(0);
}

void exithelp()
{
  COUT <<  std::endl;
  COUT << " Syntax:"<< std::endl;
  COUT <<  std::endl;
  COUT << "      eventcombiner [-v] [-i] [-n number] [-u] [-h]  inputfile1 inputfile2 ..."<< std::endl;
  COUT << " e.g  eventcombiner -v /export/rcfdata/dcm_data/built_evt/rc_3612.prdfz /export/rcfdata/dcm_data/rc/*3612*" << std::endl;
  COUT << " will combine all the *3612* (from one run number) together in the file. "<< std::endl;
  COUT << " Options:" << std::endl;
  COUT << "  -v   verbose, without that it does its work silently" << std::endl;
  COUT << "  -i   identify the events as they are processed, good for debugging" << std::endl;
  COUT << "  -e <event number> start with event number (number in evt. header)" << std::endl;
  COUT << "  -c <event number> start with nth event" << std::endl;
  COUT << "  -n <number>   stop after so many events" << std::endl;
  COUT << "  -u  write uncompressed data, default is compressed "<< std::endl;
  COUT << "  -f  force output file overwrite, normally you cannot overwrite an existing file (safety belt)"<< std::endl;
  COUT << "  -x  ignore event numbers (allow non-matching evt nrs to be combined, DANGEROUS)"<< std::endl;
  COUT << "  -h  this message" << std::endl;
  exit(0);
}

void evtcountexitmsg()
{
  COUT << "** cannot specify both -e and -c!" << std::endl;
  COUT << "    type  eventcombiner -h  for more help" << std::endl;
  exit(0);
}

oBuffer *ob;
int fd;
int file_open = 0;

int 
main(int argc, char *argv[])
{
  int c;
  int status;

  int eventnumber =0;
  int countnumber =0;
  int forceflag =0;
  int verbose = 0;
  int identify = 0;
  int maxevents = 0;
  int eventnr = 0;
  int gzipcompress = 1;
  int ignoreeventnr = 0;
  extern char *optarg;
  extern int optind;

  PHDWORD  *buffer;

  // initialize the pointers to 0;
  fd = 0;
  ob = 0;

  //  COUT << "parsing input" << std::endl;
  int buffer_size = 256*1024*16;  // makes 16MB (specifies how many dwords, so *4)



  while ((c = getopt(argc, argv, "n:c:e:viufhx")) != EOF)
    {

      switch (c) 
	{

	case 'v':   // verbose
	  verbose = 1;
	  break;

	case 'i':   // identify
	  identify = 1;
	  break;

	case 'u':   // do not gzip-compress
	  gzipcompress = 0;
	  break;

	case 'f':   // do not gzip-compress
	  forceflag = 1;
	  break;

	case 'x':   // do not gzip-compress
	  ignoreeventnr = 1;
	  break;

	case 'e':
	  if ( !sscanf(optarg, "%d", &eventnumber) ) exitmsg();
	  break;
	  
	case 'c':
	  if ( !sscanf(optarg, "%d", &countnumber) ) exitmsg();
	  break;

	case 'n':   // number of events
	  if ( !sscanf(optarg, "%d", &maxevents) ) exitmsg();
	  break;

	case 'h':
	  exithelp();
	  break;

	default:
	  break;
						
	}
    }

  if (argc < 3) exitmsg();

  if ( eventnumber && countnumber) evtcountexitmsg();

  Eventiterator *it[50];
  int no_it = 0;

  int index;
  
  char filename[256];

  strcpy(filename, argv[optind]);

  // try if the output file exists
  
  fd =  open(filename, O_RDONLY | O_LARGEFILE);
  if (fd > 0 ) 
    {
      if ( ! forceflag) 
	{
	  COUT << "file " << filename << " exists - I won't override it" << std::endl;
	  COUT << " use -f to force" << std::endl;
	  exit(1);
	}
      close (fd);
    }
  for ( index =  optind+1; index < argc; index++ ) 
    {
      COUT << "reading from file " << argv[index] << std::endl;
      it[no_it++] = new fileEventiterator(argv[index], status);
      if (status) 
	{ 
	  COUT << "could not open " <<  argv[index] << std::endl;
	  exit(1);
	}


    }

    buffer = new PHDWORD [buffer_size];

    Event *evt[50];

    int go_on = 1;

    unlink (filename);
    fd = open(filename, O_WRONLY | O_CREAT | O_EXCL | O_LARGEFILE , 
		  S_IRWXU | S_IROTH | S_IRGRP );


    if ( fd < 0 ) 
      {
	COUT << "Could not open file: " <<  filename << std::endl;
	exit (1);
      }
    if (verbose) COUT << "Opened file: " <<  filename << std::endl;
    
    if (gzipcompress) 
      {
	ob = new ogzBuffer (fd, buffer, buffer_size);
      }
    else
      {
	ob = new oBuffer (fd, buffer, buffer_size);
      }
    
    int i;
    int count = 0;

    while ( ( maxevents == 0 || eventnr < maxevents)  && go_on)
      {
	int total_length = 0;
	int enr = 0;

	for (i = 0; i< no_it; i++)
	  {
	    evt[i] = it[i]->getNextEvent();
	    if ( !evt[i] ) go_on = 0;
	    else
	      {
		//evt[i]->identify();
		total_length += evt[i]->getEvtLength();
	      }
	  }
	if (! go_on ) break;

	int *out = new int[total_length];

	int nwout;
	int current = 0;
	int take_this = 1;
	
	if ( eventnumber && evt[0]->getEvtSequence() < eventnumber)
	  take_this = 0;

	if ( countnumber && count+1 < countnumber)
	  take_this = 0;
	
	for (i = 0; i< no_it; i++)
	  {
	    if (i ==0) 
	      {	
		enr = evt[i]->getEvtSequence();
		evt[i]->Copy ( out , total_length , &nwout);
		current  = nwout;
		delete evt[i];
	      }
	    else
	      {
		if (take_this ==0  || (ignoreeventnr ==0 && evt[i]->getEvtSequence() != enr ))
		  {
		    take_this = 0;
		  }
		else
		  {
		    evt[i]->Copy (  &out[current] , total_length-current , &nwout, "DATA");
		    current += nwout;
		    out[0] +=  nwout;
		    delete evt[i];
		  }
	      }
	  }
    
	if (take_this) 
	  {
	    Event *E = new A_Event(out);
	    if (identify) E->identify();
	    
	    ob->addEvent(E);
	    delete E;
	    eventnr++;
	  }
	count++;
	delete [] out;
	
	
	
      }
    
    delete ob;
    close (fd);
}
