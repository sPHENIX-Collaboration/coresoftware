#include "SingleInttInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <ffarawobjects/InttRawHitContainerv1.h>
#include <ffarawobjects/InttRawHitv1.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <set>

SingleInttInput::SingleInttInput(const std::string &name)
  : SingleStreamingInput(name)
{
  plist = new Packet *[10];
}

SingleInttInput::~SingleInttInput()
{
  delete[] plist;
}

void SingleInttInput::FillPool(const unsigned int /*nbclks*/)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    OpenNextFile();
  }
//  std::set<uint64_t> saved_beamclocks;
   while (GetSomeMoreEvents())
  {
    Event *evt = GetEventiterator()->getNextEvent();
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt = GetEventiterator()->getNextEvent();
    }
    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
    }
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 10);

    if (npackets == 10)
    {
      exit(1);
    }
    for (int i = 0; i < npackets; i++)
    {
      if (Verbosity() > 1)
      {
        plist[i]->identify();
      }
      int num_hits = plist[i]->iValue(0, "NR_HITS");
      if (Verbosity() > 1)
      {
	std::cout << "Number of Hits: " << num_hits << " for packet "
		  << plist[i]->getIdentifier() << std::endl;
      }
      std::set<uint64_t> bclk_set;
      for (int j = 0; j < num_hits; j++)
      {
	InttRawHit *newhit = new InttRawHitv1();
	int FEE = plist[i]->iValue(j, "FEE");
	uint64_t gtm_bco = plist[i]->lValue(j, "BCO");
        newhit->set_packetid(plist[i]->getIdentifier());
	newhit->set_fee(FEE);
	newhit->set_bco(gtm_bco);
	newhit->set_adc(plist[i]->iValue(j,"ADC"));
	newhit->set_amplitude(plist[i]->iValue(j,"AMPLITUDE"));
	newhit->set_chip_id(plist[i]->iValue(j,"CHIP_ID"));
	newhit->set_channel_id(plist[i]->iValue(j,"CHANNEL_ID"));
	newhit->set_word(plist[i]->iValue(j,"DATAWORD"));
	newhit->set_FPHX_BCO(plist[i]->iValue(j,"FPHX_BCO"));
	newhit->set_full_FPHX(plist[i]->iValue(j,"FULL_FPHX"));
	newhit->set_full_ROC(plist[i]->iValue(j,"FULL_ROC"));

	gtm_bco += m_Rollover[FEE];
	bclk_set.insert(gtm_bco);
	if (gtm_bco < m_PreviousClock[FEE])
	{
	  m_Rollover[FEE] += 0x10000000000;
	  gtm_bco += 0x10000000000;  // rollover makes sure our bclks are ascending even if we roll over the 40 bit counter
	}
	m_PreviousClock[FEE] = gtm_bco;
	m_BeamClockFEE[gtm_bco].insert(FEE);
	m_FEEBclkMap[FEE] = gtm_bco;
	if (Verbosity() > 2)
	{
	  std::cout << "evtno: " << EventSequence
		    << ", hits: " << j
		    << ", nr_hits: " << num_hits
		    << ", bco: 0x" << std::hex << gtm_bco << std::dec
		    << ", FEE: " << FEE << std::endl;
	}
//          plist[i]->convert();
	if (InputManager())
	{ InputManager()->AddInttRawHit(gtm_bco, newhit); }
	m_InttRawHitMap[gtm_bco].push_back(newhit);
	m_BclkStack.insert(gtm_bco);
      }
      delete plist[i];
    }
    delete evt;
  }
//  } while (m_InttRawHitMap.size() < 10 || CheckPoolDepth(m_InttRawHitMap.begin()->first));
}

void SingleInttInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "FEE")
  {
    for (const auto &bcliter : m_BeamClockFEE)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "FEM: " << feeiter << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "FEEBCLK")
  {
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << "FEE" << bcliter.first << " bclk: 0x"
                << std::hex << bcliter.second << std::dec << std::endl;
    }
  }
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_InttRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "fee: " << feeiter->get_fee()
                  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

void SingleInttInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_InttRawHitMap)
  {
    if (iter.first <= bclk)
    {
      for (auto pktiter : iter.second)
      {
        delete pktiter;
      }
      toclearbclk.push_back(iter.first);
    }
    else
    {
      break;
    }
  }
  // for (auto iter :  m_BeamClockFEE)
  // {
  //   iter.second.clear();
  // }

  for (auto iter : toclearbclk)
  {
  m_BclkStack.erase(iter);
  m_BeamClockFEE.erase(iter);
    m_InttRawHitMap.erase(iter);
  }
}

bool SingleInttInput::CheckPoolDepth(const uint64_t bclk)
{
  // if (m_FEEBclkMap.size() < 10)
  // {
  //   std::cout << "not all FEEs in map: " << m_FEEBclkMap.size() << std::endl;
  //   return true;
  // }
  for (auto iter : m_FEEBclkMap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "my bclk 0x" << std::hex << iter.second
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    if (iter.second < bclk)
    {
      if (Verbosity() > 1)
      {
        std::cout << "FEE " << iter.first << " beamclock 0x" << std::hex << iter.second
                  << " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      }
      return false;
    }
  }
  return true;
}

void SingleInttInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
//  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  return;
}

bool SingleInttInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (CheckPoolDepth(m_InttRawHitMap.begin()->first))
  {
    if (m_InttRawHitMap.size() >= 10)
    {
      return false;
    }
  }
  return true;
}

void SingleInttInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (! dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "INTT"));
if (!detNode)
{
  detNode = new PHCompositeNode("INTT");
  dstNode->addNode(detNode);
}
  InttRawHitContainer *intthitcont = findNode::getClass<InttRawHitContainer>(detNode,"INTTRAWHIT");
  if (!intthitcont)
  {
    intthitcont = new InttRawHitContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(intthitcont, "INTTRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }

}
