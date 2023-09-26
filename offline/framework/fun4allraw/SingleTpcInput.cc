#include "SingleTpcInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <ffarawobjects/TpcRawHitContainerv1.h>
#include <ffarawobjects/TpcRawHitv1.h>

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

SingleTpcInput::SingleTpcInput(const std::string &name)
  : SingleStreamingInput(name)
{
  plist = new Packet *[10];
}

SingleTpcInput::~SingleTpcInput()
{
  delete[] plist;
}

void SingleTpcInput::FillPool(const unsigned int /*nbclks*/)
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
     int m_nWaveormInFrame = plist[i]->iValue(0, "NR_WF");
      uint64_t m_nTaggerInFrame = plist[i]->lValue(0, "N_TAGGER");
//      int m_maxFEECount = plist[i]->iValue(0, "MAX_FEECOUNT");
      std::set<uint64_t> bclk_set;
uint64_t gtm_bco = std::numeric_limits<uint64_t>::max();
    for (uint64_t t = 0; t < m_nTaggerInFrame; t++)
    {
      std::cout << "bco: 0x" << std::hex << plist[i]->lValue(t, "BCO")
		<< std::dec << std::endl;
      gtm_bco = plist[i]->lValue(t, "BCO");
      std::cout << "last bco: 0x" << std::hex << plist[i]->lValue(t, "LAST_BCO")
		<< std::dec << std::endl;
	gtm_bco += m_Rollover[0];
	bclk_set.insert(gtm_bco);
	if (gtm_bco < m_PreviousClock[0])
	{
	  m_Rollover[0] += 0x10000000000;
	  gtm_bco += m_Rollover[0];  // rollover makes sure our bclks are ascending even if we roll over the 40 bit counter
	}
/*
       m_tagger_type = (uint16_t) (p->lValue(t, "TAGGER_TYPE"));
       m_is_endat = (uint8_t) (p->lValue(t, "IS_ENDAT"));
       m_is_lvl1 = (uint8_t) (p->lValue(t, "IS_LEVEL1_TRIGGER"));
       m_bco = (uint64_t) (p->lValue(t, "BCO"));
       m_lvl1_count = (uint32_t) (p->lValue(t, "LEVEL1_COUNT"));
       m_endat_count = (uint32_t) (p->lValue(t, "ENDAT_COUNT"));
       m_last_bco = (uint64_t) (p->lValue(t, "LAST_BCO"));
       m_modebits = (uint8_t) (p->lValue(t, "MODEBITS"));
*/
    }
    for (int wf = 0; wf < m_nWaveormInFrame; wf++)
    {
      TpcRawHit *newhit = new TpcRawHitv1();
      int FEE = plist[i]->iValue(wf, "FEE");
      newhit->set_bco(plist[i]->iValue(wf, "BCO"));
      newhit->set_packetid(plist[i]->getIdentifier());
      newhit->set_samples(plist[i]->iValue(wf, "SAMPLES"));
      newhit->set_fee(FEE);
      newhit->set_channel(plist[i]->iValue(wf, "CHANNEL"));
      newhit->set_sampaaddress(plist[i]->iValue(wf, "SAMPAADDRESS"));
      newhit->set_sampachannel(plist[i]->iValue(wf, "CHANNEL"));
      m_PreviousClock[FEE] = gtm_bco;
      m_BeamClockFEE[gtm_bco].insert(FEE);
      m_FEEBclkMap[FEE] = gtm_bco;
      if (Verbosity() > 2)
      {
	std::cout << "evtno: " << EventSequence
		  << ", hits: " << wf
		  << ", num waveforms: " << m_nWaveormInFrame
		  << ", bco: 0x" << std::hex << gtm_bco << std::dec
		  << ", FEE: " << FEE << std::endl;
      }
//          plist[i]->convert();
      if (InputManager())
      {
	InputManager()->AddTpcRawHit(gtm_bco, newhit);
      }
      if (m_TpcRawHitMap.find(gtm_bco) == m_TpcRawHitMap.end())
      {
	std::vector<TpcRawHit *> intthitvector;
	m_TpcRawHitMap[gtm_bco] = intthitvector;
      }
      m_TpcRawHitMap[gtm_bco].push_back(newhit);
      m_BclkStack.insert(gtm_bco);
    }
    delete plist[i];
    }
    delete evt;
  }
//  } while (m_TpcRawHitMap.size() < 10 || CheckPoolDepth(m_TpcRawHitMap.begin()->first));
}

void SingleTpcInput::Print(const std::string &what) const
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
    for (const auto &bcliter : m_TpcRawHitMap)
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

void SingleTpcInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_TpcRawHitMap)
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
    m_TpcRawHitMap.erase(iter);
  }
}

bool SingleTpcInput::CheckPoolDepth(const uint64_t bclk)
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

void SingleTpcInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
//  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  return;
}

bool SingleTpcInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (CheckPoolDepth(m_TpcRawHitMap.begin()->first))
  {
    if (m_TpcRawHitMap.size() >= 10)
    {
      return false;
    }
  }
  return true;
}

void SingleTpcInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (! dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "TPC"));
if (!detNode)
{
  detNode = new PHCompositeNode("INTT");
  dstNode->addNode(detNode);
}
  TpcRawHitContainer *tpchitcont = findNode::getClass<TpcRawHitContainer>(detNode,"TPCRAWHIT");
  if (!tpchitcont)
  {
    tpchitcont = new TpcRawHitContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(tpchitcont, "TPCRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }

}
