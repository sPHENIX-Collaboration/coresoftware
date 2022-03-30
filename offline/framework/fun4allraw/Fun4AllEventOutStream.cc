#include "Fun4AllEventOutStream.h"

#include <fun4all/Fun4AllBase.h>                // for Fun4AllBase

#include <phool/phool.h>

#include <Event/A_Event.h>
#include <Event/Event.h>
#include <Event/oEvent.h>
#include <Event/packet.h>
#include <Event/phenixTypes.h>  // for PHDWORD

#include <algorithm>  // for copy, copy_backward, max
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ost...
#include <queue>
#include <utility>  // for swap, pair

Fun4AllEventOutStream::Fun4AllEventOutStream(const std::string &name)
  : Fun4AllBase(name)
{
}

Fun4AllEventOutStream::~Fun4AllEventOutStream()
{
  delete[] evtbuf;
  delete[] plist;
  return;
}

int Fun4AllEventOutStream::resize_evtbuf(const unsigned int newsize)
{
  delete[] evtbuf;
  evtbuf_size = newsize;
  evtbuf = new PHDWORD[evtbuf_size];
  for (unsigned int i = 0; i < evtbuf_size; i++)
  {
    evtbuf[i] = 0;
  }
  return 0;
}

int Fun4AllEventOutStream::WriteEvent(Event *evt)
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
    plist = new Packet *[max_npackets];
  }
  while ((npackets = evt->getPacketList(plist, max_npackets)) >= max_npackets)
  {
    for (int i = 0; i < npackets; i++)
    {
      delete plist[i];
    }
    delete[] plist;
    //      std::cout << "max_npackets " << max_npackets << " too small, take times 2" << std::endl;
    max_npackets *= 2;
    plist = new Packet *[max_npackets];
  }
  std::map<int, boost::numeric::interval<int> >::const_iterator dropiter;
  int dropIt;
  for (int i = 0; i < npackets; i++)
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
        std::cout << "Fun4AllEventOutStream: dropping packet " << i
             << " in list with id " << plist[i]->getIdentifier() << std::endl;
      }
    }
  }
  size += 100;  // add some size for the event header
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

int Fun4AllEventOutStream::AddPacket(const int ipkt)
{
  AddPacketRange(ipkt, ipkt);
  return 0;
}

int Fun4AllEventOutStream::DropPacket(const int ipkt)
{
  DropPacketRange(ipkt, ipkt);
  return 0;
}

int Fun4AllEventOutStream::AddPacketRange(const int minpacket, const int maxpacket)
{
  add_or_remove = 1;
  boost::numeric::interval<int> newinterval;
  try
  {
    newinterval.assign(minpacket, maxpacket);
  }
  catch (std::exception &e)
  {
    std::cout << "Exception thrown: " << e.what() << std::endl;
    std::cout << "for interval[" << minpacket << "," << maxpacket << "]" << std::endl;
    std::cout << "exiting" << std::endl;
    exit(1);
  }
  addpktrange[minpacket] = newinterval;
  if (!boost::numeric::in(minpacket, newinterval))
  {
    std::cout << PHWHERE << " boost interval does not cover minpacket " << minpacket << std::endl;
    std::cout << "that is seriously wrong, exiting" << std::endl;
    exit(1);
  }
  if (!boost::numeric::in(maxpacket, newinterval))
  {
    std::cout << PHWHERE << " boost interval does not cover maxpacket " << maxpacket << std::endl;
    std::cout << "that is seriously wrong, exiting" << std::endl;
    exit(1);
  }
  return 0;
}

int Fun4AllEventOutStream::DropPacketRange(const int minpacket, const int maxpacket)
{
  add_or_remove = 1;
  default_addall = 1;
  boost::numeric::interval<int> newinterval;
  try
  {
    newinterval.assign(minpacket, maxpacket);
  }
  catch (std::exception &e)
  {
    std::cout << "Exception thrown: " << e.what() << std::endl;
    std::cout << "for interval[" << minpacket << "," << maxpacket << "]" << std::endl;
    std::cout << "exiting" << std::endl;
    exit(1);
  }
  droppktrange[minpacket] = newinterval;
  if (!boost::numeric::in(minpacket, newinterval))
  {
    std::cout << PHWHERE << " boost interval does not cover minpacket " << minpacket << std::endl;
    std::cout << "that is seriously wrong, exiting" << std::endl;
    exit(1);
  }
  if (!boost::numeric::in(maxpacket, newinterval))
  {
    std::cout << PHWHERE << " boost interval does not cover maxpacket " << maxpacket << std::endl;
    std::cout << "that is seriously wrong, exiting" << std::endl;
    exit(1);
  }
  return 0;
}
