#include "InttCombinedRawDataDecoder.h"
#include "InttMap.h"

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
  , m_calibinfoDAC({"INTT_DACMAP", CDB})
  , m_calibinfoBCO({"INTT_BCOMAP", CDB})
{
}

int InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cerr << PHWHERE << "\n"
              << "\tCould not retrieve topNode; doing nothing" << std::endl;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator dst_itr(topNode);
  PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cerr << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)"
              << "\tCould not retrieve dst_node; doing nothing" << std::endl;
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
      std::cout << PHWHERE << "\n"
                << "\tMaking TrkrHitSetContainer" << std::endl;
    }

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
	  if (Verbosity())
      {
        std::cout << PHWHERE << "\n"
                  << "\tMaking node INTT" << std::endl;
      }

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


  ///////////////////////////////////////
  if(Verbosity())
  {
    std::cout<<"calibinfo DAC : "<<m_calibinfoDAC.first<<" "<<(m_calibinfoDAC.second==CDB?"CDB":"FILE")<<std::endl;
  }
  m_dacmap.Verbosity(Verbosity());
  if(m_calibinfoDAC.second == CDB){
     m_dacmap.LoadFromCDB(m_calibinfoDAC.first);
  } else {
     m_dacmap.LoadFromFile(m_calibinfoDAC.first);
  }
  
  ///////////////////////////////////////
  if(Verbosity())
  {
    std::cout<<"calibinfo BCO : "<<m_calibinfoBCO.first<<" "<<(m_calibinfoBCO.second==CDB?"CDB":"FILE")<<std::endl;
  }
  m_bcomap.Verbosity(Verbosity());
  if(m_calibinfoBCO.second == CDB){
     m_bcomap.LoadFromCDB(m_calibinfoBCO.first);
  } else {
     m_bcomap.LoadFromFile(m_calibinfoBCO.first);
  }

  if(m_feemap.LoadFromCDB())
  {
    std::cerr << PHWHERE << "\n"
              << "Failed to load " << m_feemap.DefaultCDBName() << " from CDB\n"
              << "Exiting" << std::endl;
	exit(1);
	gSystem->Exit(1);
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)
{
  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cerr << PHWHERE << "\n"
              << "Could not get \"TRKR_HITSET\" from Node Tree"
              << "Exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cerr << PHWHERE << "\n"
              << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << "\n"
              << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }

  if(!inttcont->get_nhits()) // empty event
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  uint64_t gl1bco = 0U;
  if(!m_runStandAlone)
  {
    Gl1RawHit* gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
	if(!gl1)
	{
      std::cerr << PHWHERE << "\n"
                << "\tCould not get node \"\", but not running standalone\n"
				<< "\tExiting" << std::endl;
	  exit(1);
	  gSystem->Exit(1);
	}
	gl1bco = gl1->get_bco();
  }
  else
  {
    gl1bco = inttcont->get_hit(0)->get_bco();
  }
  gl1bco <<= 24U;
  gl1bco >>= 24U;

  if (m_writeInttEventHeader)
  {
    intt_event_header = findNode::getClass<InttEventInfo>(topNode, "INTTEVENTHEADER");
    assert(intt_event_header);
    intt_event_header->set_bco_full(gl1bco);
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

	if(m_badmap.IsBad(raw))
	{
      continue;
	}

    int adc = intthit->get_adc();
    // amp = intthit->get_amplitude();
    uint64_t bco_full = intthit->get_bco();
    int      bco      = intthit->get_FPHX_BCO();
    
    ////////////////////////
    // bco filter
    if (m_bcomap.IsBad(raw, bco_full, bco))
    {
      continue;
    }

	if(m_feemap.Convert(ofl, raw))
	{
      std::cerr << PHWHERE << "\n"
	            << "\tconversion failed with\n"
                << "\t" << raw << "\n"
				<< "\tContinuing" << std::endl;
	  continue;
	}
    // hit_key = InttDefs::genHitKey(ofl.strip_y, ofl.strip_x);  // col, row <trackbase/InttDefs.h>
    hit_key = InttDefs::genHitKey(ofl.strip_z, ofl.strip_phi);  // col, row <trackbase/InttDefs.h>
    int time_bucket = intthit->get_bco() - gl1bco;
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

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::LoadHotChannelMapLocal(std::string const& filename)
{
  return m_badmap.LoadFromFile(filename);
}

int InttCombinedRawDataDecoder::LoadHotChannelMapRemote(std::string const& name)
{
  return m_badmap.LoadFromCDB(name);
}

