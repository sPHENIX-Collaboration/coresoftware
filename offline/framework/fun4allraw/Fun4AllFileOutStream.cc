#include "Fun4AllFileOutStream.h"

#include <fun4all/Fun4AllServer.h>

#include <Event/Event.h>
#include <Event/oBuffer.h>  // for oBuffer
#include <Event/olzoBuffer.h>

#include <phool/phool.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>  // for close
#include <cstdio>    // for snprintf
#include <cstdlib>   // for exit
#include <cstring>
#include <iostream>

Fun4AllFileOutStream::Fun4AllFileOutStream(const std::string &frule, const std::string &name)
  : Fun4AllEventOutStream(name)
  , m_FileRule(frule)
{
  memset(m_xb, 0, sizeof(m_xb));
}

Fun4AllFileOutStream::~Fun4AllFileOutStream()
{
  delete m_ob;
  if (m_OutFileDesc >= 0)
  {
    close(m_OutFileDesc);
  }
  return;
}

int Fun4AllFileOutStream::WriteEventOut(Event *evt)
{
  if (!m_ob)
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    int irun = evt->getRunNumber();
    unsigned filenamesize = m_FileRule.size() + 15;  // %010d-%04d is 14 + /0 = 15

    char *outfilename = new char[filenamesize];
    m_iSeq = se->SegmentNumber();
    int snprintfbytes = snprintf(outfilename, filenamesize, m_FileRule.c_str(), irun, m_iSeq);
    if (static_cast<unsigned>(snprintfbytes) > filenamesize)
    {
      std::cout << PHWHERE << " " << Name() << ": filename exceeds length " << filenamesize
                << ", tried " << snprintfbytes
                << ". probably it is the filerule" << m_OutFileDesc
                << " which uses other than %010d-%04d for runnumber/segment" << std::endl;
      exit(1);
    }
    m_OutFileDesc = open(outfilename, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                         S_IRWXU | S_IROTH | S_IRGRP);
    if (m_OutFileDesc == -1)  // failure to open
    {
      std::cout << "could not open " << outfilename << " quitting" << std::endl;
      exit(1);
    }
    std::cout << "opening new file " << outfilename << std::endl;
    m_ob = new olzoBuffer(m_OutFileDesc, m_xb, LENGTH, irun, m_iSeq);
    delete[] outfilename;
  }

  int status = m_ob->addEvent(evt);
  if (status)
  {
    std::cout << Name() << ": ERROR WRITING OUT FILTERED EVENT "
              << evt->getEvtSequence() << " FOR RUN "
              << evt->getRunNumber() << " Status: " << status << std::endl;
  }
  //  m_BytesWritten += 4*evt->getEvtLength(); // evtlength is in 32bit words
  m_BytesWritten = m_ob->getBytesWritten();
  if (m_BytesWritten >= m_MaxSize)
  {
    DeleteoBuffer();
    m_iSeq++;
    m_BytesWritten = 0;
    close(m_OutFileDesc);
    m_OutFileDesc = -1;
  }
  return 0;
}

int Fun4AllFileOutStream::CloseOutStream()
{
  DeleteoBuffer();
  return 0;
}

void Fun4AllFileOutStream::identify(std::ostream &os) const
{
  os << "Fun4AllFileOutStream writing to " << m_OutFileDesc << std::endl;
  return;
}

void Fun4AllFileOutStream::DeleteoBuffer()
{
  delete m_ob;
  m_ob = nullptr;
}
