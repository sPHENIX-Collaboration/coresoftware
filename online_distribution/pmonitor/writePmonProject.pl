#! /usr/bin/perl

if (@ARGV < 1) {
   print "Usage: writePmonProject.pl <projectname>\n" ;
   print "   e.g writePmonProject.pl MyAnalysis \n" ;
   exit 1;
}

$projectname=$ARGV[0];

print "creating project $projectname\n";

open (CC, "> $projectname.cc");

print CC <<EOF;

#include <iostream>
#include <pmonitor/pmonitor.h>
#include "$projectname.h"

#include <TH1.h>
#include <TH2.h>

int init_done = 0;

using namespace std;

//TH1F *h1; 
//TH2F *h2; 


int pinit()
{

  if (init_done) return 1;
  init_done = 1;

  // h1 = new TH1F ( "h1","test histogram", 100, -50.5, 49.5); 
  // h2 = new TH2F ( "h2","test histogram 2D", 100, -50.5, 49.5, 100, -500, 500);  

  return 0;

}

int process_event (Event * e)
{

  Packet *p = e->getPacket(1003);
  if (p)
    {

      //  h1->Fill ( p->iValue(0) );
      //  h2->Fill ( p->iValue(0), p->iValue(1) );

      delete p;

    }
  return 0;
}

EOF

close CC;

open (HH, "> $projectname.h");


print HH "#ifndef __";
print HH uc($projectname);
print HH "_H__\n";
print HH "#define __";
print HH uc($projectname);
print HH "_H__\n";

print HH <<EOF;

#include <Event/Event.h>

int process_event (Event *e); //++CINT 

EOF
print HH "#endif /* __";
print HH uc($projectname);
print HH "_H__ */\n";


close HH;

open (MF, "> $projectname.Makefile");

print MF <<EOF;
PACKAGE = $projectname

ROOTFLAGS = \$(shell root-config --cflags)
ROOTLIBS = \$(shell root-config --glibs)


CXXFLAGS = -I.  \$(ROOTFLAGS) -I\$(ONLINE_MAIN)/include -I\$(OFFLINE_MAIN)/include
RCFLAGS = -I.  -I\$(ONLINE_MAIN)/include -I\$(OFFLINE_MAIN)/include

LDFLAGS = -Wl,--no-as-needed  -L\$(ONLINE_MAIN)/lib -L\$(OFFLINE_MAIN)/lib -lpmonitor -lEvent -lNoRootEvent -lmessage  \$(ROOTLIBS) -fPIC



HDRFILES = \$(PACKAGE).h
LINKFILE = \$(PACKAGE)LinkDef.h


ADDITIONAL_SOURCES = 
ADDITIONAL_LIBS = 


SO = lib\$(PACKAGE).so

\$(SO) : \$(PACKAGE).cc \$(PACKAGE)_dict.C \$(ADDITIONAL_SOURCES) \$(LINKFILE)
	\$(CXX) \$(CXXFLAGS) -o \$@ -shared  \$<  \$(ADDITIONAL_SOURCES) \$(PACKAGE)_dict.C \$(LDFLAGS)  \$(ADDITIONAL_LIBS)


\$(PACKAGE)_dict.C : \$(HDRFILES) \$(LINKFILE)
	rootcint -f \$@  -c \$(RCFLAGS) \$^


.PHONY: clean

clean: 
	rm -f \$(SO) \$(PACKAGE)_dict.C \$(PACKAGE)_dict.h

EOF

close MF;

$lfname = $projectname . "LinkDef.h";
open (LF, "> $lfname");
print LF <<EOF;
#ifdef __CINT__

#pragma link C++ defined_in "$projectname.h";

#endif /* __CINT__ */
EOF

close LF;





