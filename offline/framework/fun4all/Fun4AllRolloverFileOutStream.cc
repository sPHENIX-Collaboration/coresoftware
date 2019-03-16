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
  m_Offset = offset;
  m_CurrentSequence = offset;
   m_MaxFileFize = sizeInMB;
   m_MaxFileFize =  m_MaxFileFize * 1024 * 1024;
  if (  m_MaxFileFize == 0 ||  m_MaxFileFize > MaxSize())
    {
      if (  m_MaxFileFize > MaxSize())
        {
          unsigned long long maxmb = MaxSize() / (1024 * 1024);
          cout << "setting maximum size to current max (in MB): " << maxmb << endl;
        }
       m_MaxFileFize = MaxSize();
    }
  m_Increment = increment;
  if ( m_Increment <= 0)
    {
      m_Increment = 1;  //safety belt against overwriting files
    }
}

int
Fun4AllRolloverFileOutStream::WriteEventOut(Event *evt)
{
  if (! GetoBuffer())
    {
      int irun = evt->getRunNumber();
      unsigned filenamesize = FileRule().size() + 15;

      char *outfilename = new char[filenamesize];
      iSeq(m_CurrentSequence);
      int snprintfbytes = snprintf(outfilename, filenamesize, FileRule().c_str(), irun, iSeq());
      if (static_cast<unsigned>(snprintfbytes) > filenamesize)
	{
	  cout << PHWHERE << " " << Name() << ": filename exceeds length " << filenamesize
	       << ", tried " << snprintfbytes
	       << ". probably it is the filerule" << FileRule() 
	       << " which uses other than %010d-%04d for runnumber/segment" << endl;
	  exit(1);
	}
      m_CurrentSequence += m_Increment;
      OutFileDescriptor(open(outfilename, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE ,
			     S_IRWXU | S_IROTH | S_IRGRP ));
      if (OutFileDescriptor() == -1) // failure to open
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
      SetoBuffer(new ogzBuffer ( OutFileDescriptor(), xb(), LENGTH, 6, irun, iSeq()));
      delete [] outfilename;
    }

  int status = GetoBuffer()->addEvent(evt);
  if (status)
    {
      cout << Name() << ": ERROR WRITING OUT FILTERED EVENT "
	   << evt->getEvtSequence() << " FOR RUN "
	   << evt->getRunNumber() << " Status: " << status << endl;
    }
  BytesWritten(GetoBuffer()->getBytesWritten());
  if (BytesWritten() >=  m_MaxFileFize)
    {
      DeleteoBuffer();
      BytesWritten(0);
      close(OutFileDescriptor());
      OutFileDescriptor(-1);
    }
  return 0;
}

void
Fun4AllRolloverFileOutStream::identify(ostream &os) const
{
  os << "Fun4AllRolloverFileOutStream writing to " << FileRule()
     << " current sequence " << m_CurrentSequence << endl;
    return ;
}
