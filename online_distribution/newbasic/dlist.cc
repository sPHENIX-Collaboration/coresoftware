

#include <stdlib.h>
#include <signal.h>
#include "fileEventiterator.h"
#include "testEventiterator.h"

#ifndef WIN32
#include "oncsEventiterator.h"
#endif

#include <stdio.h>

#ifdef HAVE_GETOPT_H
#ifndef WIN32
#include "getopt.h"
#else
#include "win_getopt.h"
#endif
#endif

#define ETEVENTITERATOR 1
#define FILEEVENTITERATOR 2
#define TESTEVENTITERATOR 3
#define ONCSEVENTITERATOR 4

#if defined(SunOS) || defined(Linux) || defined(OSF1)
void sig_handler(int);
#else
void sig_handler(...);
#endif

// we make it a global variable so the signal; handler can get at it.

Eventiterator *it;


void exitmsg()
{
  COUT << "** usage: dlist -ecntfTOiIh datastream" << std::endl;
  COUT << "    type  dlist -h  for more help" << std::endl;
  exit(0);
}

void evtcountexitmsg()
{
  COUT << "** cannot specify both -e and -c!" << std::endl;
  COUT << "    type  dlist -h  for more help" << std::endl;
  exit(0);
}



void exithelp()
{
  COUT << std::endl;
  int ii = 77;
  COUT << "test output " << 5 <<  ii << std::endl;

 COUT << "  dlist lists the packets contained in a given event. The event can come" << std::endl;
  COUT << "  from any of the standard sources, ET pool, file, or test stream. " << std::endl;
  COUT << "  The default is to get the next available event from a ET pool." << std::endl;
  COUT << "  dlist ONLINE would read the next event from a ET pool called ONLINE." << std::endl;
  COUT << "  dlist can optionally identify the event with the  -i option." << std::endl;
  COUT << "  You can request a certain event number with the -e option. " << std::endl;
  COUT << "  You can ask for a certain event type (data, beg-run, etc) with -t." << std::endl;
  COUT << "    -t takes a number (e.g. 9 for begin-run) for an exact match," << std::endl;
  COUT << "    or DATA or SPECIAL to selct data or special events (DATA is default)" << std::endl;
  COUT << "    you can abbreviate DATA to D or d and SPECIAL to S or s." << std::endl;
  COUT << "  Example:" << std::endl;
  COUT <<  std::endl;
  COUT << "  -T says the input is a test stream (which always has 3 packets):" << std::endl;
  COUT <<  std::endl;
  COUT << "     dlist -T" << std::endl;
  COUT << "     Packet  1001    26  0 (Unformatted)  30006 (ID4EVT)" << std::endl;
  COUT << "     Packet  1002    16  0 (Unformatted)  30005 (ID2EVT)" << std::endl;
  COUT << "     Packet  1003    10  0 (Unformatted)  30006 (ID4EVT)" << std::endl;
  COUT << std::endl;
  COUT << "  The -i option causes the event to be identified. " << std::endl;
  COUT <<  std::endl;
  COUT << "     dlist -T -i" << std::endl;
  COUT << "     -- Event     1 Run:  1331 length:    68 frames:   1 type:  1 (Data Event)" << std::endl;
  COUT << "     Packet  1001    26  0 (Unformatted)  30006 (ID4EVT)" << std::endl;
  COUT << "     Packet  1002    16  0 (Unformatted)  30005 (ID2EVT)" << std::endl;
  COUT << "     Packet  1003    10  0 (Unformatted)  30006 (ID4EVT)" << std::endl;
  COUT << std::endl;
  COUT << "  List of options: " << std::endl;
  COUT << "  usage: dlist -ecntfTOih datastream" << std::endl;
  COUT << "  -e <event number>" << std::endl;
  COUT << "  -c <number> get nth event (-e gives event with number n)" << std::endl;
  COUT << "  -n <number> repeat for n events (0: until end of stream)" << std::endl;
  COUT << "  -t <event type>" << std::endl;
  COUT << "  -i <print event identity>" << std::endl;
  COUT << "  -I <print in-depth packet identity>" << std::endl;
  COUT << "  -f (stream is a file)" << std::endl;
  COUT << "  -T (stream is a test stream)" << std::endl;
  COUT << "  -O (stream is a legacy old-style ONCS format file)" << std::endl;
  COUT << "  -v verbose" << std::endl;
  COUT << "  -h this message" << std::endl << std::endl;
  exit(0);
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

int 
main(int argc, char *argv[])
{
  int c;
  int status;

  int eventnumber =0;
  int countnumber =0;
  int repeatcount =1;
  int subeventid =0;

  int eventtypemin = 1;  
  int eventtypemax = 7;  
  int eventtype = 1;
  int type_is_a_range = 1;

  int verbose = 0;

  int ittype = ETEVENTITERATOR;
  int identify = 0;
  int fullidentify = 0;
  int generic = 0;

  extern char *optarg;
  extern int optind;

  if (argc < 2) exitmsg();

#ifndef WIN32
  while ((c = getopt(argc, argv, "n:c:e:t:iIfhTOv")) != EOF)
    switch (c) 
      {
      case 'e':
	if ( !sscanf(optarg, "%d", &eventnumber) ) exitmsg();
	break;

      case 'c':
	if ( !sscanf(optarg, "%d", &countnumber) ) exitmsg();
	break;

      case 'n':
	if ( !sscanf(optarg, "%d", &repeatcount) ) exitmsg();
	break;

      case 't':
	if (*optarg == 'S' || *optarg == 's' )
	  {
	    type_is_a_range =1;
	    eventtypemin=8;
	    eventtypemax=20;
	  }
	else if (*optarg == 'D' || *optarg == 'd')
	  {
	    type_is_a_range =1;
	    eventtypemin=1;
	    eventtypemax=7;
	  }
	else if (*optarg == 'A' || *optarg == 'a')
	  {
	    type_is_a_range =1;
	    eventtypemin=1;
	    eventtypemax=20;
	  }
	else
	  {  
	    if ( !sscanf(optarg, "%d", &eventtype) ) exitmsg();
	    type_is_a_range =0;
	  }
	break;

      case 'i':
	identify = 1;
	break;

      case 'I':
	fullidentify = 1;
	break;

       case 'T':
	ittype = TESTEVENTITERATOR;
	break;

      case 'f':
	ittype = FILEEVENTITERATOR;
	break;

      case 'O':
	ittype = ONCSEVENTITERATOR;
	break;

      case 'v':   // verbose
	verbose++;
	break;

      case 'h':
	exithelp();
	break;
      }

#else
  char* pszParam;             // gotten parameter
  char chOpt;

  while ( (chOpt = GetOption(argc, argv, "n:c:e:t:iIfhT", &pszParam) ) >1)
    {
      // COUT << "option is " << chOpt  << std::endl;
      // chOpt is valid argument
      switch (chOpt)
	{

      case 'e':
	if ( !sscanf(pszParam, "%d", &eventnumber) ) exitmsg();
	break;

      case 'c':
	if ( !sscanf(pszParam, "%d", &countnumber) ) exitmsg();
	break;

      case 'n':
	if ( !sscanf(pszParam, "%d", &repeatcount) ) exitmsg();
	break;

      case 't':
	if ( !sscanf(pszParam, "%d", &eventtype) ) exitmsg();
	break;

      case 'i':
	identify = 1;
	break;

      case 'I':
	fullidentify = 1;
	break;

       case 'T':
	ittype = TESTEVENTITERATOR;
	break;

      case 'f':
	ittype = FILEEVENTITERATOR;
	break;

      case 'h':
	exithelp();
	break;

	   
	 }
     }
#endif


  if ( eventnumber && countnumber) evtcountexitmsg();

#ifndef WIN32

  signal(SIGKILL, sig_handler);
  signal(SIGTERM, sig_handler);
  signal(SIGINT,  sig_handler);

#endif

  // see if we can open the file
  it = 0;

  switch (ittype)
    {
    case  TESTEVENTITERATOR:
      
      it = new testEventiterator();
      status =0;
      break;

    case  FILEEVENTITERATOR:
#ifndef WIN32
      it = new fileEventiterator(argv[optind], status);
#else
      //      COUT <<  " filename  is " << pszParam << std::endl;      
      it = new fileEventiterator(pszParam, status);
#endif
      break;
       
#ifndef WIN32
    case  ONCSEVENTITERATOR:
      it = new oncsEventiterator(argv[optind], status);
      break;
#else
      status = 1;
      break;
#endif
    default:
      exitmsg();
      break;
    }

  if (status)
    {
      delete it;
      it = 0;
      COUT << "Could not open input stream" << std::endl;
      exit(1);
    }

  // set the verbosity of the iterator
  it->setVerbosity(verbose);
  

  // ok. now go through the events
  Event *evt;
  Packet *s;
  int run;

  evt = it->getNextEvent();

  // first check if our requested event number is already past:

  //  if ( eventnumber && evt->getEvtSequence() > eventnumber)
  //  {
  //    COUT << "requested Event number " << eventnumber<< 
  //	" but found " <<  evt->getEvtSequence() << " already" << std::endl;
  //   delete evt;
  //   evt = 0;
  // }

  int take_this;
  int i;

  Packet *p[10000];
  int count = 0;
  int accepted_count = 0;

  while (evt)    // as long as there is a next event...
    {
      take_this = 1;
      count++;

      if ( eventnumber )
	{
	  if ( evt->getEvtSequence() == eventnumber)
	    eventnumber = 0;
	  else
	    take_this = 0;
	}

      if ( countnumber && count < countnumber)
	take_this = 0;


      if ( type_is_a_range)
	{
	  if ( (evt->getEvtType() < eventtypemin) || 
	       (evt->getEvtType() > eventtypemax) 
	       ) take_this = 0;
	}
      else
	{
	  if (evt->getEvtType() != eventtype) take_this = 0;
	}

      if (take_this)
	{
	  if (identify) evt->identify();

	  int nw = evt->getPacketList(p, 10000);
	  for (i=0; i<nw; i++)
	    {
	      if (fullidentify) p[i]->fullIdentify();
	      else p[i]->identify();
	      delete p[i];
	    }
	  delete evt;
	  if ( repeatcount==0 ||  ++accepted_count < repeatcount) evt = it->getNextEvent();
	  else	evt = 0;
	}
      else
	{
	  delete evt;
	  evt = it->getNextEvent(); // try to get the next event
	}
    }

  delete it;
  it = 0;

  return 0;
}


