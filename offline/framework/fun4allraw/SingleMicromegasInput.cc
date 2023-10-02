#include "SingleMicromegasInput.h"

#include "Fun4AllEvtInputPoolManager.h"

#include <ffarawobjects/MicromegasRawHitContainerv1.h>
#include <ffarawobjects/MicromegasRawHitv1.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <memory>
#include <set>

//______________________________________________________________
SingleMicromegasInput::SingleMicromegasInput(const std::string &name)
: SingleStreamingInput(name)
{
  plist = new Packet *[10];
}

//______________________________________________________________
SingleMicromegasInput::~SingleMicromegasInput()
{ delete[] plist; }

//______________________________________________________________
void SingleMicromegasInput::FillPool(const unsigned int /*nbclks*/)
{
  
  std::cout << "SingleMicromegasInput::FillPool" << std::endl;
  
  if (AllDone())  // no more files and all events read
  { return; }
  
  while( !GetEventiterator() )  // at startup this is a null pointer
  { OpenNextFile(); }
  
  while (GetSomeMoreEvents())
  {
    std::unique_ptr<Event> evt( GetEventiterator()->getNextEvent() );
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      
      // get next event
      evt.reset( GetEventiterator()->getNextEvent() );
    }
    
    if (Verbosity() > 2)
    { std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl; }
    
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    { evt->identify(); }
    
    if (evt->getEvtType() != DATAEVENT)
    { m_NumSpecialEvents++; }
    
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 10);
    
    if (npackets == 10)
    { exit(1); }
    
    for (int i = 0; i < npackets; i++)
    {
      
      // keep pointer to local packet
      auto& packet = plist[i];
      
      // get packet id
      const auto packet_id = packet->getIdentifier();
        
      if (Verbosity() > 1)
      { packet->identify(); }
      
      int m_nWaveormInFrame = packet->iValue(0, "NR_WF");

      // by default use previous bco clock for gtm bco
      auto& previous_bco = m_packet_bco[packet_id];
      uint64_t gtm_bco = previous_bco;

      // loop over taggers
      const int ntagger = packet->lValue(0, "N_TAGGER");
      for (int t = 0; t < ntagger; t++)
      {
        
        // only store gtm_bco for level1 type of taggers (not ENDDAT)
        const auto is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
        if( is_lvl1 )
        {
          
          gtm_bco = packet->lValue(t, "BCO");
          // std::cout << "SingleMicromegasInput::FillPool - bco: 0x" << std::hex << gtm_bco << std::dec << std::endl;
          std::cout << "SingleMicromegasInput::FillPool - bco: " << gtm_bco << std::endl;
          
          // store
          previous_bco = gtm_bco;
          
        }
      }
      
      for (int wf = 0; wf < m_nWaveormInFrame; wf++)
      {
        
        auto newhit = new MicromegasRawHitv1();
        int FEE = packet->iValue(wf, "FEE");
        newhit->set_bco(packet->iValue(wf, "BCO"));
        
        // store gtm bco in hit
        newhit->set_gtm_bco(gtm_bco);

        // packet id, fee id, channel, etc.
        newhit->set_packetid(packet_id);
        newhit->set_fee(FEE);
        newhit->set_channel(packet->iValue(wf, "CHANNEL"));
        newhit->set_sampaaddress(packet->iValue(wf, "SAMPAADDRESS"));
        newhit->set_sampachannel(packet->iValue(wf, "CHANNEL"));

//         // checksum and checksum error
//         newhit->set_checksum( packet->iValue(iwf, "CHECKSUM") );
//         newhit->set_checksum_error( packet->iValue(iwf, "CHECKSUMERROR") );

        // samples
        const uint16_t samples = packet->iValue(wf, "SAMPLES");
        newhit->set_samples( samples );

        // adc values
        for( uint16_t is =0; is < samples; ++is )
        { newhit->set_adc( is, packet->iValue( wf, is ) ); }
        
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
          
        if (InputManager())
        { InputManager()->AddMicromegasRawHit(gtm_bco, newhit); }
        
        m_MicromegasRawHitMap[gtm_bco].push_back(newhit);
        m_BclkStack.insert(gtm_bco);
      }
      
      delete packet;
    }
  }
  
  std::cout << "SingleMicromegasInput::FillPool - done." << std::endl;

}

//______________________________________________________________
void SingleMicromegasInput::Print(const std::string &what) const
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
    for (const auto &bcliter : m_MicromegasRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout
          << "fee: " << feeiter->get_fee() 
          << " at " << std::hex << feeiter << std::dec
          << std::endl;
      }
    }
  }
  
  if (what == "ALL" || what == "STACK")
  {
    for( const auto& iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}
    
//____________________________________________________________________________
void SingleMicromegasInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_MicromegasRawHitMap)
  {
    if (iter.first <= bclk)
    {
      for (auto pktiter : iter.second)
      { delete pktiter; }
      
      toclearbclk.push_back(iter.first);
    } else {
      break;
    }
    
  }
  
  for (auto iter : toclearbclk)
  {
    m_BclkStack.erase(iter);
    m_BeamClockFEE.erase(iter);
    m_MicromegasRawHitMap.erase(iter);
  }
}

//____________________________________________________________________________
bool SingleMicromegasInput::CheckPoolDepth(const uint64_t bclk)
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
      std::cout 
        << "my bclk 0x" << std::hex << iter.second
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

//_______________________________________________________
void SingleMicromegasInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
  //  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  return;
}

//_______________________________________________________
bool SingleMicromegasInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (CheckPoolDepth(m_MicromegasRawHitMap.begin()->first))
  {
    if (m_MicromegasRawHitMap.size() >= 10)
    {
      return false;
    }
  }
  return true;
}

//_______________________________________________________
void SingleMicromegasInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (! dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  
  PHNodeIterator iterDst(dstNode);
  auto detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "MICROMEGAS"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("MICROMEGAS");
    dstNode->addNode(detNode);
  }
  
  auto container = findNode::getClass<MicromegasRawHitContainer>(detNode,"MICROMEGASRAWHIT");
  if (!container)
  {
    container = new MicromegasRawHitContainerv1();
    auto newNode = new PHIODataNode<PHObject>(container, "MICROMEGASRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }
  
}
