#include "InttCombinedRawDataDecoder.h"

#include <trackbase/InttDefs.h>
#include <trackbase/InttEventInfov1.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitsetkey
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/Gl1Packet.h>
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

#include <cstdlib>     // for exit
#include <filesystem>  // for filesystem::exist
#include <format>
#include <iostream>    // for operator<<, endl, bas...
#include <map>         // for _Rb_tree_iterator

InttCombinedRawDataDecoder::InttCombinedRawDataDecoder(std::string const& name)
  : SubsysReco(name)
  , m_calibinfoDAC({"INTT_DACMAP", CDB})
  , m_calibinfoBCO({"INTT_BCOMAP", CDB})
{
  // Do nothing
  // Consider calling LoadHotChannelMapRemote()
}

int InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cout
      << PHWHERE "\n"
      << "\tCould not retrieve topNode; doing nothing\n"
	  << std::flush;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator dst_itr(topNode);
  PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cout
      << PHWHERE << "\n"
      << "\tCould not retrieve dst_node; doing nothing\n"
      << std::flush;
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
    if (Verbosity())
    {
      std::cout
        << PHWHERE << "\n"
        << "\tMaking TrkrHitSetContainer\n"
        << std::flush;
    }

    trkr_hit_set_container = new TrkrHitSetContainerv1;
    PHIODataNode<PHObject>* new_node = new PHIODataNode<PHObject>(trkr_hit_set_container, "TRKR_HITSET", "PHObject");
    trkr_node->addNode(new_node);
  }

  // Check if INTT event header already exists
  if (m_writeInttEventHeader)
  {
    auto *inttNode = dynamic_cast<PHCompositeNode*>(trkr_itr.findFirst("PHCompositeNode", "INTT"));
    if (!inttNode)
    {
      inttNode = new PHCompositeNode("INTT");
      dst_node->addNode(inttNode);
    }

    intt_event_header = findNode::getClass<InttEventInfo>(inttNode, "INTTEVENTHEADER");
    if (!intt_event_header)
    {
      intt_event_header = new InttEventInfov1();
      auto *newHeader = new PHIODataNode<PHObject>(intt_event_header, "INTTEVENTHEADER", "PHObject");
      inttNode->addNode(newHeader);
    }
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout
      << PHWHERE << "\n"
      << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree\n"
      << "removing module\n"
      << std::flush;

    Fun4AllServer* se = Fun4AllServer::instance();
    se->unregisterSubsystem(this);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  ///////////////////////////////////////
  std::cout << "calibinfo DAC : " << m_calibinfoDAC.first << " " << (m_calibinfoDAC.second == CDB ? "CDB" : "FILE") << std::endl;
  m_dacmap.Verbosity(Verbosity());
  if (m_calibinfoDAC.second == CDB)
  {
    m_dacmap.LoadFromCDB(m_calibinfoDAC.first);
  }
  else
  {
    m_dacmap.LoadFromFile(m_calibinfoDAC.first);
  }

  ///////////////////////////////////////
  std::cout << "calibinfo BCO : " << m_calibinfoBCO.first << " " << (m_calibinfoBCO.second == CDB ? "CDB" : "FILE") << std::endl;
  m_bcomap.Verbosity(Verbosity());
  int temp_offset = 0;
  if (m_calibinfoBCO.second == CDB)
  {
    temp_offset = m_bcomap.LoadFromCDB(m_calibinfoBCO.first);
  }
  else
  {
    temp_offset = m_bcomap.LoadFromFile(m_calibinfoBCO.first);
  }
  if(m_triggeredMode)
  {
    set_inttFeeOffset(temp_offset);
  }

  /// If user hasn't called with custom calibration, load default
  if (!m_badmap.OfflineLoaded() && !m_badmap.RawDataLoaded())
  {
    m_badmap.Load(); // Method loads with default tag
  }
  if (Verbosity())
  {
    std::cout << "InttBadChannelMap size: " << m_badmap.size() << std::endl;
  }
  if (1 < Verbosity())
  {
    m_badmap.Print();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)
{
  evt_inttHits_vec.clear();
  evt_ChipHit_count_map.clear();

  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cout
      << PHWHERE << "\n"
      << "\tCould not get \"TRKR_HITSET\" from Node Tree\n"
      << "\tExiting\n"
      << std::flush;
    gSystem->Exit(1);
    exit(1);

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }
  // Gl1RawHit* gl1 = nullptr;
  Gl1Packet* gl1 = nullptr;
  uint64_t gl1rawhitbco = 0;
  if (!m_runStandAlone)
  {
    gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (gl1)
    {
      gl1rawhitbco = gl1->lValue(0, "BCO");
    }
    else
    {
      auto *oldgl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
      if(!oldgl1)
      {
        std::cout << PHWHERE << " no gl1 container, exiting" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      gl1rawhitbco = oldgl1->get_bco();
    }
  }

  // get the last 40 bits by bit shifting left then right to match
  // to the mvtx bco
  auto lbshift = gl1rawhitbco << 24U;  // clang-tidy: mark as unsigned
  auto gl1bco = lbshift >> 24U;        // clang-tidy: mark as unsigned

  if (m_writeInttEventHeader)
  {
    intt_event_header = findNode::getClass<InttEventInfo>(topNode, "INTTEVENTHEADER");
    assert(intt_event_header);
    if (inttcont->get_nhits() > 0)
    {
      intt_event_header->set_bco_full(inttcont->get_hit(0)->get_bco());
    }
    else
    {
      intt_event_header->set_bco_full(0);
    }
  }

  TrkrDefs::hitsetkey hit_set_key = 0;
  TrkrDefs::hitkey hit_key = 0;
  TrkrHitSetContainer::Iterator hit_set_container_itr;
  TrkrHit* hit = nullptr;

  unsigned int loop_start = (m_SaturatedChipRejection) ? 0 : 1; // note : if m_SaturatedChipRejection, then first loop for counting, otherwise, start with the second loop
  for (unsigned int loop = loop_start; loop < 2; loop++){
    
    for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
    {
      InttRawHit* intthit = inttcont->get_hit(i);
      InttNameSpace::RawData_s raw = InttNameSpace::RawFromHit(intthit);

      int adc = intthit->get_adc();
      // amp = intthit->get_amplitude();
      uint64_t bco_full = intthit->get_bco();
      int bco = intthit->get_FPHX_BCO();

      InttNameSpace::Offline_s ofl = InttNameSpace::ToOffline(raw);

      ////////////////////////
      // bad channel filter
      if (m_badmap.OfflineLoaded() && m_badmap.IsBad(ofl))
      {
        if (1 < Verbosity())
        {
          std::cout
            << PHWHERE << "\n"
            << "\tMasking channel:\n"
            << "\t" << ofl.layer << " " << ofl.ladder_phi << " " << ofl.ladder_z << " " << ofl.strip_y << " " << ofl.strip_x << "\n"
            << std::endl;
        }
        continue;
      }

      if (m_badmap.RawDataLoaded() && m_badmap.IsBad(raw))
      {
        if (1 < Verbosity())
        {
          std::cout
            << PHWHERE << "\n"
            << "\tMasking (raw) channel:\n"
            << "\t" << raw.felix_server << " " << raw.felix_channel << " " << raw.chip << " " << raw.channel << "\n"
            << std::endl;
        }
        continue;
      }

      ////////////////////////
      // bco filter
      if (m_bcomap.IsBad(raw, bco_full, bco) && m_bcoFilter)
      {
        // std::cout<<"bad bco removed : "<<raw.felix_server<<" "<<raw.felix_channel<<" "<<raw.chip<<" "<<raw.channel<<std::endl;
        continue;
      }

      int time_bucket = 0;

      // Note: Case 1: Local, Triggered, With BCO filter
      if ( m_runStandAlone && m_triggeredMode && m_bcoFilter )
      { // NOLINT (bugprone-branch-clone)
	time_bucket = 0;
      }

      // Note: Case 2: Local, Triggered, No BCO filter
      else if ( m_runStandAlone && m_triggeredMode && !m_bcoFilter )
      { // NOLINT (bugprone-branch-clone)
	time_bucket = (intthit->get_FPHX_BCO() - (intthit->get_bco() & 0x7fU) - m_inttFeeOffset + 128) % 128;
      }

      // Note: Case 3: Local, Streaming, With BCO filter
      else if ( m_runStandAlone && !m_triggeredMode && m_bcoFilter )
      {
          std::cout<< PHWHERE << "\n" << "You selected INTT local mode, streaming mode, and BCO_filter, which is not supported. Exiting."<< std::endl;
          gSystem->Exit(1);
          exit(1);
      }

      // Note: Case 4: Local, Streaming, No BCO filter
      else if ( m_runStandAlone && !m_triggeredMode && !m_bcoFilter )
      {
          std::cout<< PHWHERE << "\n"<< "You selected INTT local mode, streaming mode, and WITHOUT BCO_filter, but GL1 is not available for time-bucket calculation. Exiting."<< std::endl;
          gSystem->Exit(1);
          exit(1);
      }

      // Note: Case 5: Global, Triggered, With BCO filter
      else if ( !m_runStandAlone && m_triggeredMode && m_bcoFilter )
      { // NOLINT (bugprone-branch-clone)
	time_bucket = 0;
      }

      // Note: Case 6: Global, Triggered, No BCO filter
      else if ( !m_runStandAlone && m_triggeredMode && !m_bcoFilter )
      { // NOLINT (bugprone-branch-clone)
	time_bucket = (intthit->get_FPHX_BCO() - (intthit->get_bco() & 0x7fU) - m_inttFeeOffset + 128) % 128;
      }

      // Note: Case 7: Global, Streaming, With BCO filter
      else if ( !m_runStandAlone && !m_triggeredMode && m_bcoFilter )
      {
          std::cout<< PHWHERE << "\n"<< "You selected INTT global mode, streaming mode, and BCO_filter, which is not supported. Exiting."<< std::endl;
          gSystem->Exit(1);
          exit(1);
      }

      // Note: Case 8: Global, Streaming, No BCO filter
      else if ( !m_runStandAlone && !m_triggeredMode && !m_bcoFilter )
      {
	time_bucket = intthit->get_FPHX_BCO() + intthit->get_bco() - gl1bco -  m_inttFeeOffset;
      }

      if(m_outputBcoDiff)
      {
        int bco_diff = 0;
        if(m_triggeredMode)
          {
            bco_diff = (intthit->get_FPHX_BCO() - (intthit->get_bco() & 0x7fU) + 128) % 128;
          }
        else
          {
            bco_diff =  intthit->get_FPHX_BCO() + intthit->get_bco() - gl1bco;
          }

        std::cout << " bco: " << " fee " << intthit->get_fee() 
            << " rawhitbco " <<  intthit->get_bco() 
            << " gl1bco " << gl1bco 
            << "  intthit->get_FPHX_BCO() " <<  intthit->get_FPHX_BCO()
            << " bcodiff " << bco_diff 
            << " time_bucket " << time_bucket 
            << std::endl;
      }

      int used_timing = (m_bcoFilter) ? intthit->get_FPHX_BCO() : time_bucket;
      std::string text_to_hit =  std::format("{}_{}_{}_{}_{}",int(intthit->get_packetid() - 3001), int(intthit->get_fee()), ((intthit->get_chip_id() - 1) % 26), int(intthit->get_channel_id()), used_timing);
      std::string text_to_chip = std::format("{}_{}_{}_{}",   int(intthit->get_packetid() - 3001), int(intthit->get_fee()), ((intthit->get_chip_id() - 1) % 26), used_timing);

      if (loop == 0) { // note : the first "loop" is for counting the number of chip hit, so we don't touch the hit set key things 

        if (std::find(evt_inttHits_vec.begin(), evt_inttHits_vec.end(),text_to_hit) == evt_inttHits_vec.end()){ // note : this is a new hit to the evt_inttHits_vec
          evt_inttHits_vec.push_back(text_to_hit);
          
          if (!evt_ChipHit_count_map.contains(text_to_chip)){
            evt_ChipHit_count_map[text_to_chip] = 1;
          }
          else{
            evt_ChipHit_count_map[text_to_chip] += 1;
          }
        }
      
        continue;
      }

      if (loop == 1 && m_SaturatedChipRejection && evt_ChipHit_count_map[text_to_chip] >= HighChipMultiplicityCut){
        
        if (Verbosity() > 10000){
          std::cout << PHWHERE << " hit in Saturated chip removed, the hit: "<<text_to_hit<<", the chip: " << text_to_chip << " with " << evt_ChipHit_count_map[text_to_chip] << " hits "<< std::endl;
        }
      
        continue;
      }

      hit_key = InttDefs::genHitKey(ofl.strip_y, ofl.strip_x);  // col, row <trackbase/InttDefs.h>

      hit_set_key = InttDefs::genHitSetKey(ofl.layer, ofl.ladder_z, ofl.ladder_phi, time_bucket);
      hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
      hit = hit_set_container_itr->second->getHit(hit_key);

      if (hit)
      {
        continue;
      }

      ////////////////////////
      // dac conversion
      int dac = m_dacmap.GetDAC(raw, adc);

      hit = new TrkrHitv2;
      //--hit->setAdc(adc);
      hit->setAdc(dac);
      hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
    }

  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}
