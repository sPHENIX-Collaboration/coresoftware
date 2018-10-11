#include "Fun4AllFileOutStream.h"
#include "Fun4AllServer.h"

#include <Event/olzoBuffer.h>
#include <Event/Event.h>

#include <phool/phool.h>

#include <cstring>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;


Fun4AllFileOutStream::Fun4AllFileOutStream(const string &frule, const string &name): 
  Fun4AllEventOutStream(name),
  filerule(frule),
  ob(NULL),
  iseq(0),
  outfile_desc(-1),
  byteswritten(0),
  MAXSIZE(10000000000LL)
{
  memset(xb, 0, sizeof(xb));
}

Fun4AllFileOutStream::~Fun4AllFileOutStream()
{
  if (ob)
    {
      delete ob;
    }
  if (outfile_desc >= 0)
    {
      close(outfile_desc);
    }
  return;
}

int
Fun4AllFileOutStream::WriteEventOut(Event *evt)
{
  if (! ob)
    {
      Fun4AllServer *se = Fun4AllServer::instance();
      int irun = evt->getRunNumber();
      unsigned filenamesize = filerule.size() + 15; // %010d-%04d is 14 + /0 = 15

      char *outfilename = new char[filenamesize];
      iseq = se->SegmentNumber();
      int snprintfbytes = snprintf(outfilename, filenamesize, filerule.c_str(), irun, iseq);
      if (static_cast<unsigned>(snprintfbytes) > filenamesize)
	{
	  cout << PHWHERE << " " << Name() << ": filename exceeds length " << filenamesize
	       << ", tried " << snprintfbytes
	       << ". probably it is the filerule" << filerule 
	       << " which uses other than %010d-%04d for runnumber/segment" << endl;
	  exit(1);
	}
      outfile_desc = open(outfilename, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE ,
                          S_IRWXU | S_IROTH | S_IRGRP );
      if (outfile_desc == -1) // failure to open
	{
	  cout << "could not open " << outfilename << " quitting" << endl;
	  exit(1);
	}
      cout << "opening new file " << outfilename << endl;
      ob = new olzoBuffer ( outfile_desc, xb, LENGTH, irun, iseq);
      delete [] outfilename;
    }

  int status = ob->addEvent(evt);
  if (status)
    {
      cout << Name() << ": ERROR WRITING OUT FILTERED EVENT "
           << evt->getEvtSequence() << " FOR RUN "
           << evt->getRunNumber() << " Status: " << status << endl;
    }
  //  byteswritten += 4*evt->getEvtLength(); // evtlength is in 32bit words
  byteswritten = ob->getBytesWritten();
  if (byteswritten >= MAXSIZE)
    {
      delete ob;
      ob = 0;
      iseq++;
      byteswritten = 0;
      close(outfile_desc);
      outfile_desc = -1;
    }
  return 0;
}

int
Fun4AllFileOutStream::CloseOutStream()
{
  if (ob)
    {
      delete ob;
      ob = 0;
    }
  return 0;
}

void 
Fun4AllFileOutStream::identify(ostream &os) const
{
  os << "Fun4AllFileOutStream writing to " << filerule << endl;
  return;
}
