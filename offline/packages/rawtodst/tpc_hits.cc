
#include "tpc_hits.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitSetv1.h>
#include <trackbase/TrkrHitv2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>


//____________________________________________________________________________..
tpc_hits::tpc_hits(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "tpc_hits::tpc_hits(const std::string &name)" << std::endl;
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;
  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv");
}

//____________________________________________________________________________..
tpc_hits::~tpc_hits()
{
  std::cout << "tpc_hits::~tpc_hits() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int tpc_hits::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "tpc_hits::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::InitRun(PHCompositeNode *topNode)
{

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "tpc_hits::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cout << "tpc_hits::InitRun - creating TRKR_HITSET." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    // create container and add to the tree
    m_hits = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(m_hits, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }  
  topNode->print();

  // we reset the BCO for the new run
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;

  //m_hits = new TrkrHitSetContainerv1();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::process_event(PHCompositeNode *topNode)
{
  std::cout << "tpc_hits::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;

  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == nullptr)
  {
    std::cout << "tpc_hits::Process_Event - Event not found" << std::endl;
    return -1;
  }
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  Packet *p = _event->getPacket(4000);
  // check all possible TPC packets that we need to analyze
  for (int packet = 4000; packet<=4230; packet+=10)
  {
    std::cout << "tpc_hits:: Packet: "<< packet << " processing" << std::endl;
    Packet *p_tmp = _event->getPacket(packet);
  
    if (!p_tmp)
    {
      std::cout << "tpc_hits:: Event getPacket: "<< packet << " IS NOT FOUND" << std::endl;
      //return Fun4AllReturnCodes::DISCARDEVENT;
    }else{
      //p = _event->getPacket(packet);
      std::cout << "tpc_hits:: Event getPacket: "<< packet << " FOUND!!!" << std::endl;
  
    }
  }
  int nr_of_waveforms = p->iValue(0, "NR_WF");

  for (auto &l : m_hitset)
  {
    l = new TrkrHitSetv1();

    int wf;
    for (wf = 0; wf < nr_of_waveforms; wf++)
    {
      int current_BCO = p->iValue(wf, "BCO") + rollover_value;

      if (starting_BCO < 0)
      {
        starting_BCO = current_BCO;
      }

      if (current_BCO < starting_BCO)  // we have a rollover
      {
        rollover_value += 0x100000;
        current_BCO = p->iValue(wf, "BCO") + rollover_value;
        starting_BCO = current_BCO;
        current_BCOBIN++;
      }
      int sampa_nr = p->iValue(wf, "SAMPAADDRESS");
      int channel = p->iValue(wf, "CHANNEL");

      int fee = p->iValue(wf, "FEE");
      int sampaAddress = p->iValue(wf, "SAMPAADDRESS");
      int sampaChannel = p->iValue(wf, "SAMPACHANNEL");
      int checksum = p->iValue(wf, "CHECKSUM");
      int checksumError = p->iValue(wf, "CHECKSUMERROR");
      int samples = p->iValue( wf, "SAMPLES" );
      


      //std::cout << "tpc_hits::Process_Event SAMPAADDRESS " << sampa_nr << std::endl;
      //std::cout << "tpc_hits::Process_Event Chn " << channel << std::endl;
      int layer;
      if (channel < 128)
      {
        layer = sampa_nr * 2;
      }
      else
      {
        layer = sampa_nr * 2 + 1;
      }
      m_hit = new TrkrHitv2();
      //      mhit->setAdc(
      int sector = 0;
      int side = 0;
      //std::cout << "tpc_hits::Process_Event Chn 1" << std::endl;
      TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(layer, sector, side);
      //std::cout << "tpc_hits::Process_Event Chn 2" << std::endl;
      TrkrHitSetContainer::Iterator hitsetit = m_hits->findOrAddHitSet(tpcHitSetKey);

      
      std::cout << "tpc_hits::Process_Event Samples "<< samples <<"Chn:"<< channel 
      <<" layer: " << layer 
      << " sampa: "<< sampa_nr 
      << " fee: "<< fee 
      << " sampaAddress: "<< sampaAddress 
      << " sampaChannel: "<< sampaChannel 
      << " checksum: "<< checksum 
      << " checksumError: "<< checksumError 

      << " hitsetkey "<< tpcHitSetKey 

      << std::endl;
      //for( int is = 0 ; is < samples ; is++ )
      //{
      //  p->iValue(wf,is);
      //}

      //std::cout << "tpc_hits::Process_Event Chn 3" << std::endl;
      for (int s = 0; s < samples; s++)
      {
        unsigned short pad = 0;
        unsigned short t = s + 2 * (current_BCO - starting_BCO);
        std::cout << "current_BCO - starting_BCO="<< current_BCO <<"-"<< starting_BCO<<"=" << s + 2 * (current_BCO - starting_BCO) << std::endl;
        // generate hit key
        //TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(pad, t);//s + 2 * (current_BCO - starting_BCO));
        TrkrDefs::hitkey hitkey = TpcDefs::genHitKey((unsigned int) pad, (unsigned int) t);
        // find existing hit, or create
        std::cout << "| " << hitkey << " "<< p->iValue(wf,s) ;
        m_hit = hitsetit->second->getHit(hitkey);
      //  // create hit, assign adc and insert in hitset     
        m_hit->setAdc(0+s);
      //  hitsetit->second->addHitSpecificKey(hitkey, hit);
      }
      std::cout << std::endl;

    }
  }
  // we skip the mapping to real pads at first. We just say
  // that we get 16 rows (segment R2) with 128 pads
  // so each FEE fills 2 rows. Not right, but one step at a time.

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::ResetEvent(PHCompositeNode * /*topNode*/)
{
  std::cout << "tpc_hits::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::EndRun(const int runnumber)
{
  std::cout << "tpc_hits::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "tpc_hits::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::Reset(PHCompositeNode * /*topNode*/)
{
  std::cout << "tpc_hits::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void tpc_hits::Print(const std::string &what) const
{
  std::cout << "tpc_hits::Print(const std::string &what) const Printing info for " << what << std::endl;
}
