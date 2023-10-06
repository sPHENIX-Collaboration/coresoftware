#include "SingleMvtxInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <ffarawobjects/MvtxRawHitContainerv1.h>
#include <ffarawobjects/MvtxRawHitv1.h>

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
//    int EventSequence = evt->getEvtSequence();
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

//          auto hbfSize = plist[i]->iValue(feeId, "NR_HBF");
          auto num_strobes = plist[i]->iValue(feeId, "NR_STROBES");
//          auto num_L1Trgs = plist[i]->iValue(feeId, "NR_PHYS_TRG");
//          for ( int iL1 = 0; iL1 < num_L1Trgs; ++iL1 )
//          {
//            os << "L1: " << iL1  << std::hex << " BCO: 0x" << lValue(feeId, iL1, "L1_IR_BCO");
//            os << std::dec << " BC: " << iValue(feeId, iL1, "L1_IR_BC") << endl;
//          }

          for (int i_strb{0}; i_strb < num_strobes; ++i_strb)
          {
            //TODO: need to add virtual function lvalue(const int, const int , char*)
            // using 0 by now
            auto strb_bco = 0; //plist[i]->lValue(feeId, i_strb, "TRG_IR_BCO");
            auto strb_bc  = plist[i]->iValue(feeId, i_strb, "TRG_IR_BC");
            auto num_hits = plist[i]->iValue(feeId, i_strb, "TRG_NR_HITS");
            //	if (Verbosity() > 2)
            //	{
            //	  std::cout << "evtno: " << EventSequence
            //		    << ", hits: " << j
            //		    << ", nr_hits: " << num_hits
            //		    << ", bco: 0x" << std::hex << gtm_bco << std::dec
            //		    << ", FEE: " << FEE << std::endl;
            //	}

            for (int i_hit{0}; i_hit < num_hits; ++i_hit)
            {
              MvtxRawHit *newhit = new MvtxRawHitv1();
//              newhit->set_bco(strb_bco);
              newhit->set_strobe_bc(strb_bc);
              newhit->set_chip_bc(plist[i]->iValue(feeId, i_strb, i_hit, "HIT_BC"));
//            TODO: get layer,stave, chipId from feeId and iValue(feeId, i_trg, i_hit, "HIT_CHIP_ID");
//              newhit->set_layer_id();
//              newhit->set_stave_id();
//              newhit->set_chip_id();
              newhit->set_row(plist[i]->iValue(feeId, i_strb, i_hit, "HIT_ROW"));
              newhit->set_col(plist[i]->iValue(feeId, i_strb, i_hit, "HIT_COL"));
              if (InputManager())
              {
                InputManager()->AddMvtxRawHit(strb_bco, newhit);
              }
              m_MvtxRawHitMap[strb_bco].push_back(newhit);
            }
            m_BclkStack.insert(strb_bco);
          }
        }
      }
      plist[i]->convert();
      delete plist[i];
    }
    delete evt;
  }
//  } while (m_MvtxRawHitMap.size() < 10 || CheckPoolDepth(m_MvtxRawHitMap.begin()->first));
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
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
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

void SingleMvtxInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
//  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  return;
}

bool SingleMvtxInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (CheckPoolDepth(m_MvtxRawHitMap.begin()->first))
  {
    if (m_MvtxRawHitMap.size() >= 10)
    {
      return false;
    }
  }
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

  MvtxRawHitContainer* mvtxhitcont = findNode::getClass<MvtxRawHitContainer>(detNode,"MVTXRAWHIT");
  if (! mvtxhitcont)
  {
    mvtxhitcont = new MvtxRawHitContainerv1();
    PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(mvtxhitcont, "MVTXRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }
}
