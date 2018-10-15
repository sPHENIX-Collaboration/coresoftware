#include "Fun4AllEventOutputManager.h"
#include "Fun4AllServer.h"
#include "Fun4AllEventOutStream.h"
#include "Fun4AllRolloverFileOutStream.h"

#include <phool/getClass.h>

#include <Event/Event.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

Fun4AllEventOutputManager::Fun4AllEventOutputManager(const string &myname
						     , const string &filerule
						     , const unsigned int sizeInMB
						     , const int offset
						     , const int increment ):
  Fun4AllOutputManager( myname )
{
  outfilerule = filerule;
  outstream = new Fun4AllRolloverFileOutStream(filerule, sizeInMB, offset, increment);
  outstream->SetManager(this);
  return ;
}

Fun4AllEventOutputManager::~Fun4AllEventOutputManager()
{
  if (outstream)
    {
      delete outstream;
    }
  return ;
}


void
Fun4AllEventOutputManager::Print(const string &what) const
{
  cout << Name() << " writes " << outfilerule << endl;
  // base class print method
  Fun4AllOutputManager::Print( what );

  return ;
}

int Fun4AllEventOutputManager::Write(PHCompositeNode* /*startNode*/)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  PHCompositeNode *topNode = se->topNode();
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (!evt)
    {
      cout << PHWHERE << "0 Event Pointer" << endl;
      return -1;
    }
  outstream->WriteEvent(evt);
  return 0;
}

int
Fun4AllEventOutputManager::AddPacket(const int ipkt)
{
  int iret = -1;
  if (outstream)
    {
      iret = outstream->AddPacket(ipkt);
    }
  else
    {
      cout << PHWHERE << "Cannot add packet" << endl;
    }
  return iret;
}

int
Fun4AllEventOutputManager::AddPacketRange(const int ipktmin, const int ipktmax)
{
  int iret = -1;
  if (outstream)
    {
      iret = outstream->AddPacketRange(ipktmin, ipktmax);
    }
  else
    {
      cout << PHWHERE << "Cannot add packet" << endl;
    }
  return iret;
}

int
Fun4AllEventOutputManager::DropPacket(const int ipkt)
{
  int iret = -1;
  if (outstream)
    {
      iret = outstream->DropPacket(ipkt);
    }
  else
    {
      cout << PHWHERE << "Cannot drop packet" << endl;
    }
  return iret;
}

int
Fun4AllEventOutputManager::DropPacketRange(const int ipktmin, const int ipktmax)
{
  int iret = -1;
  if (outstream)
    {
      iret = outstream->DropPacketRange(ipktmin, ipktmax);
    }
  else
    {
      cout << PHWHERE << "Cannot drop packet" << endl;
    }
  return iret;
}

void
Fun4AllEventOutputManager::SetOutfileName(const std::string &fname)
{
  outfilename = fname;
  return ;
}
