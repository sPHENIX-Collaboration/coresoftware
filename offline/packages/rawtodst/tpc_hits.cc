
#include "tpc_hits.h"

#include <trackbase/TpcDefs.h>
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

#include <phool/PHCompositeNode.h>

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
  // std::cout << "tpc_hits::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;

  topNode->print();

  // we reset the BCO for the new run
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;

  m_hits = new TrkrHitSetContainerv1();

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

  Packet *p = _event->getPacket(4001);
  if (!p)
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
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
      TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layer, sector, side);
      TrkrHitSetContainer::Iterator hitsetit = m_hits->findOrAddHitSet(hitsetkey);
      for (int s = 0; s < p->iValue(wf, "SAMPLES"); s++)
      {
        int pad = 0;
        TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(pad, s + 2 * (current_BCO - starting_BCO));
        TrkrHit *hit = hitsetit->second->getHit(hitkey);
        hit->setAdc(0);
      }
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
