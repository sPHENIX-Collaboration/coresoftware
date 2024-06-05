#include "InttCombinedRawDataDecoder.h"
#include "InttMapping.h"

#include <trackbase/InttDefs.h>
#include <trackbase/InttEventInfov1.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitsetkey
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cstdlib>   // for exit
#include <filesystem>// for filesystem::exist
#include <iostream>  // for operator<<, endl, bas...
#include <map>       // for _Rb_tree_iterator

InttCombinedRawDataDecoder::InttCombinedRawDataDecoder(std::string const& name)
  : SubsysReco(name)
{
  // Do nothing
  // Consider calling LoadHotChannelMapRemote()
}

int InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\tCould not retrieve top node; doing nothing" << std::endl;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator dst_itr(topNode);
  PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\tCould not retrieve dst node; doing nothing" << std::endl;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator trkr_itr(dst_node);
  PHCompositeNode* trkr_node = dynamic_cast<PHCompositeNode*>(trkr_itr.findFirst("PHCompositeNode", "TRKR"));
  if (!trkr_node)
  {
    trkr_node = new PHCompositeNode("TRKR");
    dst_node->addNode(trkr_node);
  }

  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    trkr_hit_set_container = new TrkrHitSetContainerv1;
    PHIODataNode<PHObject>* new_node = new PHIODataNode<PHObject>(trkr_hit_set_container, "TRKR_HITSET", "PHObject");
    trkr_node->addNode(new_node);
  }

  // Check if INTT event header already exists
  if (m_writeInttEventHeader)
  {
    auto inttNode = dynamic_cast<PHCompositeNode *>(trkr_itr.findFirst("PHCompositeNode", "INTT"));
    if (!inttNode)
    {
      inttNode = new PHCompositeNode("INTT");
      dst_node->addNode(inttNode);
    }

    intt_event_header = findNode::getClass<InttEventInfo>(inttNode, "INTTEVENTHEADER");
    if (!intt_event_header)
    {
      intt_event_header = new InttEventInfov1();
      auto newHeader = new PHIODataNode<PHObject>(intt_event_header, "INTTEVENTHEADER", "PHObject");
      inttNode->addNode(newHeader);
    }
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "removing module" << std::endl;

    Fun4AllServer* se = Fun4AllServer::instance();
    se->unregisterSubsystem(this);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Check if member calibrations have been loaded at this point
  // (user may have already called some of the Load... at the macro level)
  // If not, try to load them (with default names)

  if(!m_badmap.IsLoaded())
  {
    if(LoadBadMap() && Verbosity())
	{
      std::cout << __PRETTY_FUNCTION__ << "\n"
                << "\tCould not load m_badmap\n"
                << "\tNo hot channel filtering will be performed\n"
                << "\tDecoder will still run" << std::endl;
	}
  }

  if(!m_bcomap.IsLoaded())
  {
    if(LoadBcoMap() && Verbosity())
	{
      std::cout << __PRETTY_FUNCTION__ << "\n"
                << "\tCould not load m_bcomap\n"
                << "\tNo time window filetering will be performed\n"
                << "\tDecoder will still run" << std::endl;
	}
  }

  if(!m_dacmap.IsLoaded())
  {
    if(LoadDacMap() && Verbosity())
	{
      std::cout << __PRETTY_FUNCTION__ << "\n"
                << "\tCould not load m_dacmap\n"
                << "\tDefault digital to analog conversion will be used\n"
                << "\tDecoder will still run" << std::endl;
	}
  }

  if(!m_feemap.IsLoaded())
  {
    if(LoadFeeMap()) // FeeMap is special; we need to load it or we cannot assign RAWHIT info to ladders
    {
      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\tUnable to load m_feemap\n"
                << "\tCannot map rawhits to ladders\n"
                << "\tExiting" << std::endl;
  	  gSystem->Exit(1);
  	  exit(1);
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)
{
  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cerr << PHWHERE << "\n"
              << "Could not get \"TRKR_HITSET\" from Node Tree\n"
              << "Exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cerr << PHWHERE << "\n"
              << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree"
              << "Exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  Gl1RawHit* gl1 = nullptr;
  if (!m_runStandAlone)
  {
    gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
    if (!gl1)
    {
      std::cout << PHWHERE << " no gl1 container, exiting" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  uint64_t gl1rawhitbco = m_runStandAlone ? 0 : gl1->get_bco();
  // get the last 40 bits by bit shifting left then right to match
  // to the mvtx bco
  auto lbshift = gl1rawhitbco << 24U; // clang-tidy: mark as unsigned
  auto gl1bco = lbshift >> 24U; // clang-tidy: mark as unsigned

  if (m_writeInttEventHeader)
  {
    intt_event_header = findNode::getClass<InttEventInfo>(topNode, "INTTEVENTHEADER");
    assert(intt_event_header);
    intt_event_header->set_bco_full(inttcont->get_hit(0)->get_bco());
  }

  TrkrDefs::hitsetkey hit_set_key = 0;
  TrkrDefs::hitkey hit_key = 0;
  TrkrHitSetContainer::Iterator hit_set_container_itr;
  TrkrHit* hit = nullptr;

  InttMap::RawData_s raw;
  InttMap::Offline_s ofl;
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit* intthit = inttcont->get_hit(i);

    raw.pid = intthit->get_packetid();
    raw.fee = intthit->get_fee();
    raw.chp = (intthit->get_chip_id() + 25) % 26;
    raw.chn = intthit->get_channel_id();

    int adc = intthit->get_adc();
    // amp = intthit->get_amplitude();
    uint64_t bco_full = intthit->get_bco();
    int      bco      = intthit->get_FPHX_BCO();


	m_feemap.Convert(ofl, raw);

	if(m_badmap.IsBad(ofl))
	{
      continue;
	}

	if(m_bcomap.IsBad(raw.pid - 3001, raw.fee, bco_full, bco))
	{
      continue;
	}

    hit_key = InttDefs::genHitKey(ofl.strip_z, ofl.strip_phi);  // col, row <trackbase/InttDefs.h>
    int time_bucket = m_runStandAlone ? 0 : intthit->get_bco() - gl1bco;
    hit_set_key = InttDefs::genHitSetKey(ofl.layer, ofl.ladder_z, ofl.ladder_phi, time_bucket);
    hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
    hit = hit_set_container_itr->second->getHit(hit_key);

    if (hit)
    {
      continue;
    }

    ////////////////////////
    // dac conversion
    int dac = m_dacmap.GetDAC(raw.pid - 3001, raw.fee, raw.chp, raw.chn, adc);

    hit = new TrkrHitv2;
    //--hit->setAdc(adc);
    hit->setAdc(dac);
    hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
