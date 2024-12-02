#include "Fun4AllRolloverFileOutStream.h"

#include "Fun4AllEventOutputManager.h"

#include <Event/Event.h>
#include <Event/oBuffer.h>  // for oBuffer
#include <Event/ophBuffer.h>

#include <phool/phool.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>  // for close
#include <cstdio>    // for snprintf
#include <cstdlib>   // for exit
#include <iostream>

Fun4AllRolloverFileOutStream::Fun4AllRolloverFileOutStream(const std::string &frule,
                                                           const unsigned int nEvents,
                                                           const unsigned int sizeInMB,
                                                           const int offset,
                                                           const int increment,
                                                           const std::string &name)
  : Fun4AllFileOutStream(frule, name)

{
  m_Offset = offset;
  m_CurrentSequence = offset;
  m_MaxNEvents = nEvents;
  m_MaxFileFize = sizeInMB;
  m_MaxFileFize = m_MaxFileFize * 1024 * 1024;
  if (m_MaxFileFize == 0 || m_MaxFileFize > MaxSize())
  {
    if (m_MaxFileFize > MaxSize())
    {
      uint64_t maxmb = MaxSize() / (1024ULL * 1024ULL);
      std::cout << "setting maximum size to current max (in MB): " << maxmb << std::endl;
    }
    m_MaxFileFize = MaxSize();
  }
  m_Increment = increment;
  if (m_Increment <= 0)
  {
    m_Increment = 1;  // safety belt against overwriting files
  }
}

int Fun4AllRolloverFileOutStream::WriteEventOut(Event *evt)
{
  if (!GetoBuffer())
  {
    int irun = evt->getRunNumber();
    unsigned filenamesize = FileRule().size() + 15;

    char *outfilename = new char[filenamesize];
    iSeq(m_CurrentSequence);
    // NOLINTNEXTLINE(hicpp-vararg)
    int snprintfbytes = snprintf(outfilename, filenamesize, FileRule().c_str(), irun, iSeq());
    if (static_cast<unsigned>(snprintfbytes) > filenamesize)
    {
      std::cout << PHWHERE << " " << Name() << ": filename exceeds length " << filenamesize
                << ", tried " << snprintfbytes
                << ". probably it is the filerule" << FileRule()
                << " which uses other than %010d-%04d for runnumber/segment" << std::endl;
      exit(1);
    }
    m_CurrentSequence += m_Increment;
    // NOLINTNEXTLINE(hicpp-vararg)
    OutFileDescriptor(open(outfilename, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                           S_IRWXU | S_IROTH | S_IRGRP));
    if (OutFileDescriptor() == -1)  // failure to open
    {
      std::cout << "could not open " << outfilename << " quitting" << std::endl;
      exit(1);
    }
    if (Verbosity() > 0)
    {
      std::cout << "Fun4AllRolloverFileOutStream: opening new file " << outfilename << std::endl;
    }
    MyManager()->SetOutfileName(outfilename);
    SetoBuffer(new ophBuffer(OutFileDescriptor(), xb(), LENGTH, irun, iSeq()));
    delete[] outfilename;
  }

  int status = GetoBuffer()->addEvent(evt);
  if (status)
  {
    std::cout << Name() << ": ERROR WRITING OUT FILTERED EVENT "
              << evt->getEvtSequence() << " FOR RUN "
              << evt->getRunNumber() << " Status: " << status << std::endl;
  }
  SetNEvents(GetNEvents() + 1);
  BytesWritten(GetoBuffer()->getBytesWritten());

  if (m_MaxNEvents > 0 && GetNEvents() >= m_MaxNEvents)
  {
    open_new_file();
  }
  if (BytesWritten() >= m_MaxFileFize)
  {
    open_new_file();
  }
  return 0;
}
void Fun4AllRolloverFileOutStream::identify(std::ostream &os) const
{
  os << "Fun4AllRolloverFileOutStream writing to " << FileRule()
     << " current sequence " << m_CurrentSequence << std::endl;
  return;
}

void Fun4AllRolloverFileOutStream::open_new_file()
{
  DeleteoBuffer();
  BytesWritten(0);
  SetNEvents(0);
  close(OutFileDescriptor());
  OutFileDescriptor(-1);
  MyManager()->RunAfterClosing();
}
