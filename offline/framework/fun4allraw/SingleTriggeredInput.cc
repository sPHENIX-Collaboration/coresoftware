#include "SingleTriggeredInput.h"

#include <frog/FROG.h>

#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <set>
#include <utility>  // for pair
#include <vector>

SingleTriggeredInput::SingleTriggeredInput(const std::string &name)
  : Fun4AllBase(name)
{
  m_bclkarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkdiffarray.fill(std::numeric_limits<uint64_t>::max());
}

SingleTriggeredInput::~SingleTriggeredInput()
{
  while (m_EventDeque.begin() != m_EventDeque.end())
  {
    delete m_EventDeque.front();
    m_EventDeque.pop_front();
  }
  delete m_EventIterator;
}

int SingleTriggeredInput::fileopen(const std::string &filenam)
{
  std::cout << PHWHERE << "trying to open " << filenam << std::endl;
  if (IsOpen())
  {
    std::cout << "Closing currently open file "
              << FileName()
              << " and opening " << filenam << std::endl;
    fileclose();
  }
  FileName(filenam);
  FROG frog;
  std::string fname = frog.location(FileName());
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << FileName() << std::endl;
  }
  int status = 0;
  m_EventIterator = new fileEventiterator(fname.c_str(), status);
  m_EventsThisFile = 0;
  if (status)
  {
    delete m_EventIterator;
    m_EventIterator = nullptr;
    std::cout << PHWHERE << Name() << ": could not open file " << fname << std::endl;
    return -1;
  }
  IsOpen(1);
  AddToFileOpened(fname);  // add file to the list of files which were opened
  return 0;
}

int SingleTriggeredInput::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  delete m_EventIterator;
  m_EventIterator = nullptr;
  IsOpen(0);
  // if we have a file list, move next entry to top of the list
  // or repeat the same entry again
  UpdateFileList();
  return 0;
}

int SingleTriggeredInput::FillEventVector()
{
  while (GetEventIterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return -1;
    }
  }
  if (!m_EventDeque.empty())
  {
    return 0;
  }
  size_t i{0};
  uint64_t tmp = m_bclkarray[pooldepth];
  m_bclkarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkdiffarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkarray[0] = tmp;
  // std::cout << std::hex << "m_bclkarray[0]: 0x" << m_bclkarray[0] << ", ind " << pooldepth
  // 	    << ", 0x" << m_bclkarray[pooldepth] << std::endl;
  while (i < pooldepth)
  {
    Event *evt = GetEventIterator()->getNextEvent();
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        FilesDone(1);
        if (Verbosity() > 0)
        {
          std::cout << "no more events to read, deque depth: " << m_EventDeque.size() << std::endl;
        }
        return -1;
      }
      evt = GetEventIterator()->getNextEvent();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      delete evt;
      continue;
    }
    //    std::cout << Name() << std::endl;
    //    evt->identify();
    evt->convert();
    if (firstcall)
    {
      CreateDSTNodes(evt);
      firstcall = false;
    }
    //    std::cout << Name() << ":filling index " << i+1 << " with " << std::hex << bla << std::dec << std::endl;
    uint64_t myClock = GetClock(evt);
    if (myClock == std::numeric_limits<uint64_t>::max())
    {
      delete evt;
      continue;
    }
    m_bclkarray[i + 1] = GetClock(evt);
    if (m_Event == 0)
    {
      m_bclkdiffarray[i] = 0;
    }
    else
    {
      if (m_bclkarray[i + 1] < m_bclkarray[i])
      {
        m_bclkdiffarray[i] = (m_bclkarray[i + 1] + 0x100000000) - m_bclkarray[i];  // handle rollover
      }
      else
      {
        m_bclkdiffarray[i] = m_bclkarray[i + 1] - m_bclkarray[i];
      }
    }
    //    std::cout << "Event " << evt->getEvtSequence() << "bclk: " << m_bclkarray[i+1] << std::endl;

    m_Event++;
    // std::cout << "m_bclkarray[" << i<< "]: 0x" << std::hex << m_bclkarray[i] << std::dec
    // 	      << ", m_bclkarray[" << i << "+1]: 0x" << std::hex << m_bclkarray[i+1] << std::dec
    // 	      << ", m_bclkdiffarray[" << i << "]: 0x" << std::hex << m_bclkdiffarray[i] << std::dec << std::endl;
    m_EventDeque.push_back(evt);

    i++;
  }
  // std::cout << Name() << std::endl;
  // for (auto iter : m_bclkdiffarray)
  // {
  //   std::cout << "0x" << std::hex << iter << std::dec << std::endl;
  // }
  return m_EventDeque.size();
}

uint64_t SingleTriggeredInput::GetClock(Event *evt)
{
  std::vector<Packet *> pktvec = evt->getPacketVector();
  uint64_t clock = static_cast<uint64_t>(pktvec[0]->lValue(0, "CLOCK") & 0xFFFFFFFF);  // NOLINT (hicpp-signed-bitwise)
  if (Verbosity() > 1)
  {
    std::cout << "Getting clock from packet " << pktvec[0]->getIdentifier() << std::endl;
  }
  //     uint64_t clock = pktvec[0]->lValue(0, "CLOCK");
  //     pktvec[0]->identify();
  // std::cout << "pkt event: " << pktvec[0]->iValue(0, "EVTNR") << ", clock: 0x"
  // 	       << std::hex << pktvec[0]->lValue(0, "CLOCK") << std::dec << std::endl;
  for (auto *iter : pktvec)
  {
    delete iter;
  }
  return clock;
}

void SingleTriggeredInput::FillPool(const unsigned int keep)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  if (!FilesDone())
  {
    if (FillEventVector() != 0)
    {
      bool isequal = std::equal(begin(), end(), Gl1Input()->begin());
      if (!isequal)
      {
        const auto *iter1 = begin();
        const auto *iter2 = Gl1Input()->begin();
        //        auto iter3 = beginclock();
        //        auto iter4 = Gl1Input()->beginclock();
        int position = 0;
        //        int ifirst = 1;
        while (iter1 != end())
        {
          // std::cout << "position " << position << " test 0x" << std::hex
          //  	      << *iter1 << " versus 0x" << *iter2 << std::dec << std::endl;
          //  std::cout  << "bclk1: 0x" << *iter3 << ", gl1bclk: 0x" << *iter4 << std::dec << std::endl;
          //	    m_EventDeque[position]->identify();
          if (*iter1 != *iter2)
          {
            if (*iter1 == std::numeric_limits<uint64_t>::max())
            {
              std::cout << Name() << " No more events found marked by uint64_t max clock" << std::endl;
              FilesDone(1);
              break;
            }
            if (Verbosity() > 0)
            {
              std::cout << "remove Event " << m_EventDeque[position]->getEvtSequence() << " clock: 0x"
                        << std::hex << GetClock(m_EventDeque[position]) << ", clkdiff(me) 0x"
                        << *iter1 << ", clkdiff(gl1) 0x" << *iter2
                        << std::dec << std::endl;
            }
            //	    std::cout << "and booom" << std::endl;
            //	    std::cout <<  "remove Event " << m_EventDeque[position]->getEvtSequence() << " clock: " << m_EventDeque[position]->getPacket(6067)->lValue(0, "CLOCK")<< std::endl ;
            // position is the first bad index,
            for (size_t i = position; i < m_EventDeque.size(); ++i)
            {
              delete m_EventDeque[i];
            }
            m_EventDeque.erase(m_EventDeque.begin() + (position), m_EventDeque.end());
            break;
          }

          // std::cout <<  "good Event " << m_EventDeque[position]->getEvtSequence() << " clock: " << m_EventDeque[position]->getPacket(6067)->lValue(0, "CLOCK")<< std::endl ;

          ++position;
          ++iter1;
          ++iter2;
          // if (ifirst)
          // {
          //   ifirst = 0;
          // }
          // else
          // {
          //   ++iter3;
          //   ++iter4;
          // }
        }
        std::cout << "Aborting event loop after processing remaining good events" << std::endl;
        FilesDone(1);
      }

      //	std::cout << "we are good" << std::endl;
    }
  }
  if (keep > 100000000)
  {
    std::cout << "huh?" << std::endl;
  }
  if (m_EventDeque.empty())
  {
    std::cout << Name() << ":all events done" << std::endl;
    AllDone(1);
    return;
  }
  if (Verbosity() > 1)
  {
    std::cout << "deque size: " << m_EventDeque.size() << std::endl;
  }
  Event *evt = m_EventDeque.front();
  m_EventDeque.pop_front();
  RunNumber(evt->getRunNumber());
  //  std::cout << "Handling event " << evt->getEvtSequence();
//  CaloPacket *newhit = new CaloPacketv1();
  std::vector<Packet *> pktvec = evt->getPacketVector();
  for (auto *packet : pktvec)
  {
    int packet_id = packet->getIdentifier();
    CaloPacket *newhit = findNode::getClass<CaloPacket>(m_topNode,packet_id);
    newhit->setPacketEvtSequence(packet->iValue(0, "EVTNR"));
    int nr_modules = packet->iValue(0, "NRMODULES");
    int nr_channels = packet->iValue(0, "CHANNELS");
    int nr_samples = packet->iValue(0, "SAMPLES");
    newhit->setNrModules(nr_modules);
    newhit->setNrChannels(nr_channels);
    newhit->setNrSamples(nr_samples);
    newhit->setIdentifier(packet_id);
    newhit->setBCO(packet->lValue(0, "CLOCK"));
    //     std::cout << ", clock :" << packet->lValue(0, "CLOCK") << std::endl;
    for (int ifem = 0; ifem < nr_modules; ifem++)
    {
      newhit->setFemClock(ifem, packet->iValue(ifem, "FEMCLOCK"));
      newhit->setFemEvtSequence(ifem, packet->iValue(ifem, "FEMEVTNR"));
      newhit->setFemSlot(ifem, packet->iValue(ifem, "FEMSLOT"));
      newhit->setChecksumLsb(ifem, packet->iValue(ifem, "CHECKSUMLSB"));
      newhit->setChecksumMsb(ifem, packet->iValue(ifem, "CHECKSUMMSB"));
      newhit->setCalcChecksumLsb(ifem, packet->iValue(ifem, "CALCCHECKSUMLSB"));
      newhit->setCalcChecksumMsb(ifem, packet->iValue(ifem, "CALCCHECKSUMMSB"));
    }
    for (int ipmt = 0; ipmt < nr_channels; ipmt++)
    {
      // store pre/post only for suppressed channels, the array in the packet routines is not
      // initialized so reading pre/post for not zero suppressed channels returns garbage
      bool isSuppressed = packet->iValue(ipmt, "SUPPRESSED");
      newhit->setSuppressed(ipmt, isSuppressed);
      if (isSuppressed)
      {
        newhit->setPre(ipmt, packet->iValue(ipmt, "PRE"));
        newhit->setPost(ipmt, packet->iValue(ipmt, "POST"));
      }
      else
      {
        for (int isamp = 0; isamp < nr_samples; isamp++)
        {
          newhit->setSample(ipmt, isamp, packet->iValue(isamp, ipmt));
        }
      }
    }
    delete packet;
  }
  delete evt;
}

void SingleTriggeredInput::CreateDSTNodes(Event *evt)
{
  static std::string CompositeNodeName = "Packets";
  PHNodeIterator iter(m_topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    m_topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", CompositeNodeName));
  if (!detNode)
  {
    detNode = new PHCompositeNode(CompositeNodeName);
    dstNode->addNode(detNode);
  }
 std::vector<Packet *> pktvec = evt->getPacketVector();
 for (auto *piter : pktvec)
  {
    int packet_id = piter->getIdentifier();
    m_PacketSet.insert(packet_id);
    std::string PacketNodeName = std::to_string(packet_id);
  CaloPacket *calopacket = findNode::getClass<CaloPacket>(detNode, PacketNodeName);
    if (!calopacket)
    {
      calopacket = new CaloPacketv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(calopacket, PacketNodeName, "PHObject");
    detNode->addNode(newNode);
  }
    delete piter;
  }
}
