#include <stdlib.h>
#include <signal.h>
#include "fileEventiterator.h"
//#include <packetConstants.h>
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
  COUT << "** usage: ddump -ecnstdfghiIFTOv datastream" << std::endl;
  COUT << "    type  ddump -h   for more help" << std::endl;
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
  COUT << "   ddump dumps selected packets from a given event. The event can come" << std::endl;
  COUT << "   from any of the standard sources, ET pool, file, or test stream. " << std::endl;
  COUT << "   The default is to get the next available event from a ET pool." << std::endl;
  COUT << std::endl;
  COUT << "   ddump can optionally identify the event with the  -i option." << std::endl;
  COUT << "   You can request a certain event number with the -e option. " << std::endl;
  COUT << "   You can ask for a certain event type (data, beg-run, etc) with -t." << std::endl;
  COUT << "     -t takes a number (e.g. 9 for begin-run) for an exact match," << std::endl;
  COUT << "     or DATA or SPECIAL to selct data or special events (DATA is default)" << std::endl;
  COUT << "     you can abbreviate DATA to D or d and SPECIAL to S or s." << std::endl;

  COUT << std::endl;
  COUT << "   Note: In order to find out which packets are contained in an event, you can" << std::endl;
  COUT << "   use the dlist utility." << std::endl;
  COUT << std::endl;
  COUT << "   Example:" << std::endl;
  COUT << std::endl;
  COUT << "   -T says that the input stream is a test stream (where events have 3 " << std::endl;
  COUT << "   predictable packets 1001, 1002, and 1003); we dump packet 1002:" << std::endl;
  COUT << std::endl;
  COUT << "   > ddump -p 1002 -T" << std::endl;
  COUT << "   Packet  1002    16  0 (Unformatted)  30005 (ID2EVT)" << std::endl;
  COUT << "       0 |     1    2    3    4    5    6    7    8 " << std::endl;
  COUT << "       8 |     9    a    b    c    d    e    f   10 " << std::endl;
  COUT << "      10 |    11   12   13   14 " << std::endl;
  COUT << std::endl;
  COUT << "   The -i option also lets the event identify itself." << std::endl;
  COUT << "   The -I option gets you an in-depth identification of the packets." << std::endl;

  COUT << std::endl;
  COUT << "   > ddump -s 1002 -T -i" << std::endl;
  COUT << "    -- Event     1 Run:  1331 length:    68 frames:   1 type:  1 (Data Event)" << std::endl;
  COUT << "   Packet  1002    16  0 (Unformatted)  30005 (ID2EVT)" << std::endl;
  COUT << "       0 |     1    2    3    4    5    6    7    8 " << std::endl;
  COUT << "       8 |     9    a    b    c    d    e    f   10 " << std::endl;
  COUT << "      10 |    11   12   13   14 " << std::endl;

  COUT << " The -F option also causes the frames in the event to be listed." << std::endl;
  COUT << " If you list a particular packet, you will get the frame where" << std::endl; 
  COUT << " this packet is in." << std::endl; 

  COUT << std::endl;

  COUT << "  List of options: " << std::endl;
  COUT << " -e <event number>" << std::endl;
  COUT << " -c <number> get nth event (-e gives event with number n)" << std::endl;
  COUT << " -n <number> repeat for n events (0: until end of stream)" << std::endl;
  COUT << " -p <Packet Id>" << std::endl;
  COUT << " -t <event type>" << std::endl;
  COUT << " -i <print event identity>" << std::endl;
  COUT << " -I <print in-depth packet identity (default is short form)>" << std::endl;
  COUT << " -f (stream is a file)" << std::endl;
  COUT << " -T (stream is a test stream)" << std::endl;
  COUT << " -O (stream is a legacy old-style ONCS format file)" << std::endl;
  COUT << " -g use generic dump" << std::endl;
  COUT << " -d numbers are std::decimal (default std::hex) for generic dump" << std::endl;
  COUT << " -o numbers are octal (default std::hex) for generic dump" << std::endl;
  COUT << " -v verbose" << std::endl;
  COUT << " -h this message" << std::endl << std::endl;
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
  int dumpstyle = EVT_HEXADECIMAL;
  int identify = 0;
  int fullidentify = 0;
  int listframe = 0;
  int listHistory = 0;
  int listError = 0;
  int generic = 0;

  extern char *optarg;
  extern int optind;

  if (argc < 2) exitmsg();

#ifndef WIN32
  while ((c = getopt(argc, argv, "n:c:e:s:p:t:idfghIFTOHEv")) != EOF)
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

      case 'p':
	if ( !sscanf(optarg, "%d", &subeventid) ) exitmsg();
	break;

      case 's':
	if ( !sscanf(optarg, "%d", &subeventid) ) exitmsg();
	break;

      case 't':
	if (*optarg == 'S' || *optarg == 's' )
	  {
	    type_is_a_range =1;
	    eventtypemin=8;
	    eventtypemax=20;
	  }
	else if (*optarg == 'A' || *optarg == 'a')
	  {
	    type_is_a_range =1;
	    eventtypemin=1;
	    eventtypemax=20;
	  }
	else if (*optarg == 'D' || *optarg == 'd')
	  {
	    type_is_a_range =1;
	    eventtypemin=1;
	    eventtypemax=7;
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

      case 'F':
	listframe = 1;
	break;

      case 'H':
	listHistory = 1;
	break;

      case 'E':
	listError = 1;
	break;

      case 'g':
	generic = 1;
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

      case 'o':
	dumpstyle = EVT_OCTAL;
	break;

      case 'd':
	dumpstyle = EVT_DECIMAL;
	break;

      case 'v':   // verbose
	verbose++;
	break;

      case 'h':
	exithelp();
	break;
      }
#else
  char* pszParam;
  char chOpt;
  while ((chOpt = GetOption(argc, argv, "n:c:e:s:p:t:idfghIFTO", &pszParam)) >1 )
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

      case 'p':
	if ( !sscanf(pszParam, "%d", &subeventid) ) exitmsg();
	break;

      case 's':
	if ( !sscanf(pszParam, "%d", &subeventid) ) exitmsg();
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

      case 'F':
	listframe = 1;
	break;

      case 'H':
	listHistory = 1;
	break;

      case 'E':
	listError = 1;
	break;

      case 'g':
	generic = 1;
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

      case 'o':
	dumpstyle = EVT_OCTAL;
	break;

      case 'd':
	dumpstyle = EVT_DECIMAL;
	break;

      case 'h':
	exithelp();
	break;
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
      it = new fileEventiterator(pszParam, status);
#endif
      break;
     
    case  ONCSEVENTITERATOR:
#ifndef WIN32
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
  int i;

  evt = it->getNextEvent();

  // first check if our requested event number is already past:

  //  if ( eventnumber && evt->getEvtSequence() > eventnumber)
  // {
  //   COUT << "requested Event number " << eventnumber<< 
  //	" but found " <<  evt->getEvtSequence() << " already" << std::endl;
  //   delete evt;
  //   evt = 0;
  // }

  int take_this;
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

	  if ( subeventid)
	    {
	      
	      if ( (s = evt->getPacket(subeventid)) )  // get the subevent
		{
		  if (listframe)  evt->listFrame (subeventid);
		  if (fullidentify) s->fullIdentify();
		  if (generic) s->gdump(dumpstyle);
		  else s->dump();
		  if (listHistory) evt->listHistory();
		  if (listError) evt->listError();
		}
	    }
	  else
	    {
	      if (listframe)  evt->listFrame();
	      Packet *p[10000];
	      int nw = evt->getPacketList(p, 10000);
	      for (i=0; i<nw; i++)
		{
		  if (fullidentify) p[i]->fullIdentify();
		  if (generic) p[i]->gdump(dumpstyle);
		  else p[i]->dump();
		  delete p[i];
		}
	      if (listHistory) evt->listHistory();
	      if (listError) evt->listError();
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

