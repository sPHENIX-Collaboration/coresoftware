#include "Fun4AllRolloverFileOutStream.h"
#include "Fun4AllEventOutputManager.h"
#include "Fun4AllServer.h"

#include <Event/ogzBuffer.h>
#include <Event/Event.h>

#include <phool/phool.h>

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

Fun4AllRolloverFileOutStream::Fun4AllRolloverFileOutStream(const string &frule,
							   const unsigned int sizeInMB,
							   const int offset,
							   const int increment,
							   const string &name) : Fun4AllFileOutStream(frule, name)

{
  i_offset = offset;
  current_sequence = offset;
  max_file_size = sizeInMB;
  max_file_size = max_file_size * 1024 * 1024;
  if ( max_file_size == 0 || max_file_size > MAXSIZE)
    {
      if ( max_file_size > MAXSIZE)
        {
          unsigned long long maxmb = MAXSIZE / (1024 * 1024);
          cout << "setting maximum size to current max (in MB): " << maxmb << endl;
        }
      max_file_size = MAXSIZE;
    }
  i_increment = increment;
  if ( i_increment <= 0)
    {
      i_increment = 1;  //safety belt against overwriting files
    }
}

int
Fun4AllRolloverFileOutStream::WriteEventOut(Event *evt)
{
  if (! ob)
    {
      int irun = evt->getRunNumber();
      unsigned filenamesize = filerule.size() + 15;

      char *outfilename = new char[filenamesize];
      iseq = current_sequence;
      int snprintfbytes = snprintf(outfilename, filenamesize, filerule.c_str(), irun, iseq);
      if (static_cast<unsigned>(snprintfbytes) > filenamesize)
	{
	  cout << PHWHERE << " " << Name() << ": filename exceeds length " << filenamesize
	       << ", tried " << snprintfbytes
	       << ". probably it is the filerule" << filerule 
	       << " which uses other than %010d-%04d for runnumber/segment" << endl;
	  exit(1);
	}
      current_sequence += i_increment;
      outfile_desc = open(outfilename, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE ,
                          S_IRWXU | S_IROTH | S_IRGRP );
      if (outfile_desc == -1) // failure to open
	{
	  cout << "could not open " << outfilename << " quitting" << endl;
	  exit(1);
	}
      if (Verbosity() > 0)
        {
          cout << "Fun4AllRolloverFileOutStream: opening new file " << outfilename << endl;
        }
      mymanager->SetOutfileName(outfilename);
      // compression level 6 is best compromize between speed and compression
      // max is 9 which is much slower but only squeezes out a few more bytes
      ob = new ogzBuffer ( outfile_desc, xb, LENGTH, 6, irun, iseq);
      delete [] outfilename;
    }

  int status = ob->addEvent(evt);
  if (status)
    {
      cout << Name() << ": ERROR WRITING OUT FILTERED EVENT "
	   << evt->getEvtSequence() << " FOR RUN "
	   << evt->getRunNumber() << " Status: " << status << endl;
    }
  byteswritten = ob->getBytesWritten();
  if (byteswritten >= max_file_size)
    {
      delete ob;
      ob = 0;
      byteswritten = 0;
      close(outfile_desc);
      outfile_desc = -1;
    }
  return 0;
}

void
Fun4AllRolloverFileOutStream::identify(ostream &os) const
{
  os << "Fun4AllRolloverFileOutStream writing to " << filerule
     << " current sequence " << current_sequence << endl;
    return ;
}
