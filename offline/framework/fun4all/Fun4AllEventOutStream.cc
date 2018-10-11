#include "Fun4AllEventOutStream.h"

#include <Event/oEvent.h>
#include <Event/A_Event.h>
#include <Event/Event.h>
#include <Event/packet.h>

#include <phool/phool.h>

#include <queue>

using namespace std;

Fun4AllEventOutStream::Fun4AllEventOutStream(const std::string &name): 
  Fun4AllBase(name),
  evtbuf(NULL),
  evtbuf_size(0),
  add_or_remove(0),
  plist(NULL),
  max_npackets(1000), // there shouldn't be more than this number of packets in a single event
  npackets(0),
  default_addall(0),
  mymanager(NULL)
{}

Fun4AllEventOutStream::~Fun4AllEventOutStream()
{
  delete [] evtbuf;
  delete [] plist;
  return ;
}

int
Fun4AllEventOutStream::resize_evtbuf(const unsigned int newsize)
{
  delete [] evtbuf;
  evtbuf_size = newsize;
  evtbuf = new PHDWORD[evtbuf_size];
  for (unsigned int i = 0;i < evtbuf_size;i++)
    {
      evtbuf[i] = 0;
    }
  return 0;
}

int
Fun4AllEventOutStream::WriteEvent(Event *evt)
{
  int iret;
  if (!add_or_remove)
    {
      iret = WriteEventOut(evt);
      return iret;
    }
  std::queue<int> savepacket;
  unsigned int size = 0;
  if (!plist)
    {
      plist = new Packet*[max_npackets];
    }
  while ((npackets = evt->getPacketList(plist, max_npackets)) >= max_npackets)
    {
      for (int i = 0; i < npackets; i++)
        {
          delete plist[i];
        }
      delete [] plist;
      //      cout << "max_npackets " << max_npackets << " too small, take times 2" << endl;
      max_npackets *= 2;
      plist = new Packet*[max_npackets];
    }
  std::map<int, boost::numeric::interval<int> >::const_iterator dropiter;
  int dropIt;
  for (int i = 0; i < npackets;i++)
    {
      int packetid = plist[i]->getIdentifier();
      if (default_addall)
        {
          dropIt = 0;
          for (dropiter = droppktrange.begin(); dropiter != droppktrange.end(); ++dropiter)
            {
              if (packetid < dropiter->first)
                {
                  // abort loop if packetid is smaller than first packet in range
                  break;
                }
              if (boost::numeric::in(packetid, dropiter->second))
                {
                  dropIt = 1;
                  break;
                }
            }
        }
      else
        {
          dropIt = 1;
        }
      for (dropiter = addpktrange.begin(); dropiter != addpktrange.end(); ++dropiter)
        {
          if (packetid < dropiter->first)
            {
              // abort loop if packetid is smaller than first packet in range
              break;
            }
          if (boost::numeric::in(packetid, dropiter->second))
            {
              dropIt = 0;
              break;
            }
        }
      if (!dropIt)
        {
          savepacket.push(i);
          size += plist[i]->getLength() + 4;
        }
      else
        {
          if (Verbosity() > 0)
            {
              cout << "Fun4AllEventOutStream: dropping packet " << i
                   << " in list with id " << plist[i]->getIdentifier() << endl;
            }
        }
    }
  size += 100; // add some size for the event header
  if (size > evtbuf_size)
    {
      // Add 10000 so we do this resize only a few times
      resize_evtbuf(size + 10000);
    }

  oEvent new_event(evtbuf, size, evt->getRunNumber(), evt->getEvtType(), evt->getEvtSequence());
  while (!savepacket.empty())
    {
      int index = savepacket.front();
      new_event.addPacket(plist[index]);
      savepacket.pop();
    }

  Event *newE = new A_Event(evtbuf);
  iret = WriteEventOut(newE);
  delete newE;
  for (int i = 0; i < npackets; i++)
    {
      delete plist[i];
    }
  return iret;
}

int
Fun4AllEventOutStream::AddPacket(const int ipkt)
{
  AddPacketRange(ipkt, ipkt);
  return 0;
}

int
Fun4AllEventOutStream::DropPacket(const int ipkt)
{
  DropPacketRange(ipkt, ipkt);
  return 0;
}

int
Fun4AllEventOutStream::AddPacketRange(const int minpacket, const int maxpacket)
{
  add_or_remove = 1;
  boost::numeric::interval<int> newinterval;
  try
    {
      newinterval.assign(minpacket, maxpacket);
    }
  catch (exception& e)
    {
      cout << "Exception thrown: " << e.what() << endl;
      cout << "for interval[" << minpacket << "," << maxpacket << "]" << endl;
      cout << "exiting" << endl;
      exit(1);
    }
  addpktrange[minpacket] = newinterval;
  if (! boost::numeric::in(minpacket, newinterval))
    {
      cout << PHWHERE << " boost interval does not cover minpacket " << minpacket << endl;
      cout << "that is seriously wrong, exiting" << endl;
      exit(1);
    }
  if (! boost::numeric::in(maxpacket, newinterval))
    {
      cout << PHWHERE << " boost interval does not cover maxpacket " << maxpacket << endl;
      cout << "that is seriously wrong, exiting" << endl;
      exit(1);
    }
  return 0;
}

int
Fun4AllEventOutStream::DropPacketRange(const int minpacket, const int maxpacket)
{
  add_or_remove = 1;
  default_addall = 1;
  boost::numeric::interval<int> newinterval;
  try
    {
      newinterval.assign(minpacket, maxpacket);
    }
  catch (exception& e)
    {
      cout << "Exception thrown: " << e.what() << endl;
      cout << "for interval[" << minpacket << "," << maxpacket << "]" << endl;
      cout << "exiting" << endl;
      exit(1);
    }
  droppktrange[minpacket] = newinterval;
  if (! boost::numeric::in(minpacket, newinterval))
    {
      cout << PHWHERE << " boost interval does not cover minpacket " << minpacket << endl;
      cout << "that is seriously wrong, exiting" << endl;
      exit(1);
    }
  if (! boost::numeric::in(maxpacket, newinterval))
    {
      cout << PHWHERE << " boost interval does not cover maxpacket " << maxpacket << endl;
      cout << "that is seriously wrong, exiting" << endl;
      exit(1);
    }
  return 0;
}

