#include "SingleMvtxInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <ffarawobjects/MvtxRawHitContainerv1.h>
#include <ffarawobjects/MvtxRawHitv1.h>
#include <ffarawobjects/MvtxRawRunHeader.h>
#include <ffarawobjects/MvtxRawEvtHeader.h>

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

SingleMvtxInput::SingleMvtxInput(const std::string &name)
  : SingleStreamingInput(name)
{
  plist = new Packet *[2];
}

SingleMvtxInput::~SingleMvtxInput()
{
  delete[] plist;
}

void SingleMvtxInput::FillPool(const unsigned int /*nbclks*/)
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
    int npackets = evt->getPacketList(plist, 2);

    if (npackets > 2)
    {
      exit(1);
    }
    for (int i = 0; i < npackets; i++)
    {
      // Ignoring packet not from MVTX detector
      if ( (plist[i]->getIdentifier() < 2001) || (plist[i]->getIdentifier() > 2052) )
      {
        continue;
      }
      if (Verbosity() > 1)
      {
        plist[i]->identify();
      }
      int num_feeId = plist[i]->iValue(-1, "NR_LINKS");
      if (Verbosity() > 1)
      {
        std::cout << "Number of feeid in RCDAQ events: " << num_feeId << " for packet "
          << plist[i]->getIdentifier() << std::endl;
      }
      if (num_feeId > 0)
      {
        for (int i_fee{0}; i_fee < num_feeId; ++i_fee)
        {
          auto feeId = plist[i]->iValue(i_fee, "FEEID");
          auto link = DecodeFeeid(feeId);
//          auto hbfSize = plist[i]->iValue(feeId, "NR_HBF");
          auto num_strobes = plist[i]->iValue(feeId, "NR_STROBES");
          auto num_L1Trgs = plist[i]->iValue(feeId, "NR_PHYS_TRG");
          for ( int iL1 = 0; iL1 < num_L1Trgs; ++iL1 )
          {
            auto l1Trg_bco = plist[i]->lValue(feeId, iL1, "L1_IR_BCO");
//            auto l1Trg_bc  = plist[i]->iValue(feeId, iL1, "L1_IR_BC");
            m_GtmL1BcoSet.emplace(l1Trg_bco);
          }

          m_FeeStrobeMap[feeId] += num_strobes;
          for (int i_strb{0}; i_strb < num_strobes; ++i_strb)
          {
            auto strb_bco = plist[i]->lValue(feeId, i_strb, "TRG_IR_BCO");
            auto strb_bc  = plist[i]->iValue(feeId, i_strb, "TRG_IR_BC");
            auto num_hits = plist[i]->iValue(feeId, i_strb, "TRG_NR_HITS");
            if (Verbosity() > 4)
            {
              std::cout << "evtno: " << EventSequence << ", Fee: " << feeId;
              std::cout << " Layer: " << link.layer << " Stave: " << link.stave;
              std::cout << " GBT: " << link.gbtid << ", bco: 0x" << std::hex << strb_bco << std::dec;
              std::cout << ", n_hits: " << num_hits << std::endl;
            }
            for (int i_hit{0}; i_hit < num_hits; ++i_hit)
            {
              auto chip_bc = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_BC");
              auto chip_id = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_CHIP_ID");
              auto chip_row = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_ROW");
              auto chip_col = plist[i]->iValue(feeId, i_strb, i_hit, "HIT_COL");
              MvtxRawHit *newhit = new MvtxRawHitv1();
              newhit->set_bco(strb_bco);
              newhit->set_strobe_bc(strb_bc);
              newhit->set_chip_bc(chip_bc);
              newhit->set_layer_id(link.layer);
              newhit->set_stave_id(link.stave);
              newhit->set_chip_id( 3 * link.gbtid + chip_id);
              newhit->set_row(chip_row);
              newhit->set_col(chip_col);
              if (InputManager())
              {
                InputManager()->AddMvtxRawHit(strb_bco, newhit);
              }
              m_MvtxRawHitMap[strb_bco].push_back(newhit);
            }
            m_BeamClockFEE[strb_bco].insert(feeId);
            m_BclkStack.insert(strb_bco);
            m_FEEBclkMap[feeId] = strb_bco;
          }
        }
      }
//      plist[i]->convert();
      delete plist[i];
    }
    delete evt;
  }
}

void SingleMvtxInput::Print(const std::string &what) const
{
  //TODO: adapt to MVTX case
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
    for (const auto &bcliter : m_MvtxRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "fee: " << feeiter->get_stave_id()
                  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "GET_NR_STROBES")
  {
    for( auto& iter : m_FeeStrobeMap )
    {
      std::cout << "Total number of strobes for feeid: " << iter.first << ", " << iter.second << std::endl;
    }
  }
}

void SingleMvtxInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_MvtxRawHitMap)
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
    m_MvtxRawHitMap.erase(iter);
  }
}

bool SingleMvtxInput::CheckPoolDepth(const uint64_t bclk)
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
      std::cout << iter.first << " my bclk 0x" << std::hex << iter.second
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    // equal case when we have more strobe with same bco
    // due not synchronization
    if (iter.second <= bclk)
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

void SingleMvtxInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
//  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  clearGtmL1BcoSet();
  return;
}

bool SingleMvtxInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
//  if (CheckPoolDepth(m_MvtxRawHitMap.begin()->first))
//  {
    if (m_MvtxRawHitMap.size() >= 200)
    {
      return false;
    }
//  }
  return true;
}

void SingleMvtxInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (! dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "MVTX"));
  if (! detNode)
  {
    detNode = new PHCompositeNode("MVTX");
    dstNode->addNode(detNode);
  }

  MvtxRawRunHeader* mvtxRH = findNode::getClass<MvtxRawRunHeader>(detNode,"MVTXRAWRUNHEADER");
  if (! mvtxRH)
  {
    mvtxRH = new MvtxRawRunHeader();
    PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(mvtxRH, "MVTXRAWRUNHEADER", "PHObject");
    detNode->addNode(newNode);
  }

  MvtxRawEvtHeader* mvtxEH = findNode::getClass<MvtxRawEvtHeader>(detNode,"MVTXRAWEVTHEADER");
  if (! mvtxEH)
  {
    mvtxEH = new MvtxRawEvtHeader();
    PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(mvtxEH, "MVTXRAWEVTHEADER", "PHObject");
    detNode->addNode(newNode);
  }

  MvtxRawHitContainer* mvtxhitcont = findNode::getClass<MvtxRawHitContainer>(detNode,"MVTXRAWHIT");
  if (! mvtxhitcont)
  {
    mvtxhitcont = new MvtxRawHitContainerv1();
    PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(mvtxhitcont, "MVTXRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }
}

void SingleMvtxInput::clearGtmL1BcoSet()
{
  m_GtmL1BcoSet.clear();
}
