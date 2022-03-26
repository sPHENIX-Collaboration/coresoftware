#include "Fun4AllEventOutputManager.h"

#include "Fun4AllEventOutStream.h"
#include "Fun4AllRolloverFileOutStream.h"

#include <fun4all/Fun4AllOutputManager.h>  // for Fun4AllOutputManager
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Event/Event.h>

#include <iostream>
#include <string>

Fun4AllEventOutputManager::Fun4AllEventOutputManager(const std::string &myname, const std::string &filerule, const unsigned int sizeInMB, const int offset, const int increment)
  : Fun4AllOutputManager(myname)
  , m_OutFileRule(filerule)
{
  m_OutStream = new Fun4AllRolloverFileOutStream(filerule, sizeInMB, offset, increment);
  m_OutStream->SetManager(this);
  return;
}

Fun4AllEventOutputManager::~Fun4AllEventOutputManager()
{
  delete m_OutStream;
  return;
}

void Fun4AllEventOutputManager::Print(const std::string &what) const
{
  std::cout << Name() << " writes " << m_OutFileRule << std::endl;
  // base class print method
  Fun4AllOutputManager::Print(what);

  return;
}

int Fun4AllEventOutputManager::Write(PHCompositeNode * /*startNode*/)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *topNode = se->topNode();
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (!evt)
  {
    std::cout << PHWHERE << "0 Event Pointer" << std::endl;
    return -1;
  }
  m_OutStream->WriteEvent(evt);
  return 0;
}

int Fun4AllEventOutputManager::AddPacket(const int ipkt)
{
  int iret = -1;
  if (m_OutStream)
  {
    iret = m_OutStream->AddPacket(ipkt);
  }
  else
  {
    std::cout << PHWHERE << "Cannot add packet" << std::endl;
  }
  return iret;
}

int Fun4AllEventOutputManager::AddPacketRange(const int ipktmin, const int ipktmax)
{
  int iret = -1;
  if (m_OutStream)
  {
    iret = m_OutStream->AddPacketRange(ipktmin, ipktmax);
  }
  else
  {
    std::cout << PHWHERE << "Cannot add packet" << std::endl;
  }
  return iret;
}

int Fun4AllEventOutputManager::DropPacket(const int ipkt)
{
  int iret = -1;
  if (m_OutStream)
  {
    iret = m_OutStream->DropPacket(ipkt);
  }
  else
  {
    std::cout << PHWHERE << "Cannot drop packet" << std::endl;
  }
  return iret;
}

int Fun4AllEventOutputManager::DropPacketRange(const int ipktmin, const int ipktmax)
{
  int iret = -1;
  if (m_OutStream)
  {
    iret = m_OutStream->DropPacketRange(ipktmin, ipktmax);
  }
  else
  {
    std::cout << PHWHERE << "Cannot drop packet" << std::endl;
  }
  return iret;
}

void Fun4AllEventOutputManager::SetOutfileName(const std::string &fname)
{
  OutFileName(fname);
  return;
}
