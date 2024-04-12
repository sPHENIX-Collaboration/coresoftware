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
{
  m_calibs = {
    {m_badmap.DefaultCDBName(), (struct calib_load_s){
      .ptr =      &m_badmap
    }},
    {m_bcomap.DefaultCDBName(), (struct calib_load_s){
      .ptr =      &m_bcomap
    }},
    {m_dacmap.DefaultCDBName(), (struct calib_load_s){
      .ptr =      &m_dacmap
    }},
    {m_feemap.DefaultCDBName(), (struct calib_load_s){
      .method =   FROM_CDB,
	  .ptr =      &m_feemap,
      .required = true
    }},
  };
}

int InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cerr << PHWHERE << "\n"
              << "\tCould not retrieve topNode\n"
              << "\tExiting" << std::endl;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator dst_itr(topNode);
  PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cerr << PHWHERE << "\n"
              << "\tCould not retrieve dst_node\n"
              << "\tExiting" << std::endl;
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
      inttNode = new PHCompositeNode("INTT");
      dst_node->addNode(inttNode);
    }

    InttEventInfo* intt_event_header = findNode::getClass<InttEventInfo>(inttNode, "INTTEVENTHEADER");
    if (!intt_event_header)
    {
      intt_event_header = new InttEventInfov1();
      auto newHeader = new PHIODataNode<PHObject>(intt_event_header, "INTTEVENTHEADER", "PHObject");
      inttNode->addNode(newHeader);
    }
  }

  if(LoadCalibs() != 0)
  {
    std::cerr << PHWHERE << "\n"
              << "\tFailed to load all calibrations\n"
              << "\tExiting" << std::endl;
	gSystem->Exit(1);
	exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)
{
  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cerr << PHWHERE << "\n"
              << "\tCould not get node \"TRKR_HITSET\"\n"
              << "\tAborting" << std::endl;
    gSystem->Exit(1);
    exit(1);
    return Fun4AllReturnCodes::ABORTRUN;
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cerr << PHWHERE << "\n"
              << "\tCould not get node \"" << m_InttRawNodeName << "\"\n"
              << "\tAborting" << std::endl;
    gSystem->Exit(1);
    exit(1);
    return Fun4AllReturnCodes::ABORTRUN;
  }

  uint64_t gl1bco = 0U;
  if (!m_runStandAlone)
  {
    Gl1RawHit* gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
    if (!gl1)
    {
      std::cerr << PHWHERE << "\n"
                << "\tCould not get node \"GL1RAWHIT\", but not running standalone\n"
                << "\tAborting\n";
      gSystem->Exit(1);
      exit(1);
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    gl1bco = gl1->get_bco();
  }
  else
  {
    gl1bco = inttcont->get_hit(0)->get_bco();
  }
  // get the last 40 bits by bit shifting left then right to match
  gl1bco <<= 24U; // clang-tidy: mark as unsigned
  gl1bco >>= 24U; // clang-tidy: mark as unsigned

  if (m_writeInttEventHeader)
  {
    InttEventInfo* intt_event_header = findNode::getClass<InttEventInfo>(topNode, "INTTEVENTHEADER");
    if (!intt_event_header)
    {
      std::cerr << PHWHERE << "\n"
                << "\tCould not get node \"INTTEVENTHEADER\", but are writing event header\n"
                << "\tAborting\n";
      gSystem->Exit(1);
      exit(1);
      return Fun4AllReturnCodes::ABORTEVENT;
    }
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

    if (m_badmap.IsBad(raw))
    {
      continue;
    }

    if (m_feemap.Convert(ofl, raw))
    {
      if(Verbosity())
      {
          std::cerr << PHWHERE << "\n"
                    << "\tconversion failed with\n"
                    << "\t" << raw << "\n"
                    << "\tContinuing" << std::endl;
      }
      continue;
    }

    int adc = intthit->get_adc();
	int dac = m_dacmap.GetDAC (
      raw.pid - 3001,
      raw.fee,
      raw.chp,
      raw.chn,
      adc
    );
	if(dac == 0xFFFFU)
	{
	  if(Verbosity())
      {
          std::cerr << PHWHERE << "\n"
                    << "\tDAC conversion lookup failed\n"
                    << "\tContinuing" << std::endl;;
      }
      continue;
    }

    // amp = intthit->get_amplitude();
    uint64_t bco_full = intthit->get_bco();
    int      bco      = intthit->get_FPHX_BCO();
	if (m_bcomap.IsBad (
      raw.pid - 3001,
	  raw.fee,
	  bco_full,
	  bco))
	{
      continue;
	}

    // hit_key = InttDefs::genHitKey(ofl.strip_y, ofl.strip_x); // col, row <trackbase/InttDefs.h>
    hit_key = InttDefs::genHitKey(ofl.strip_z, ofl.strip_phi); // col, row <trackbase/InttDefs.h>
    int time_bucket = intthit->get_bco() - gl1bco;
    hit_set_key = InttDefs::genHitSetKey(ofl.layer, ofl.ladder_z, ofl.ladder_phi, time_bucket);
    hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
    hit = hit_set_container_itr->second->getHit(hit_key);

    if (hit)
    {
      continue;
    }

    hit = new TrkrHitv2;

    hit->setAdc(dac); // adc
    hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
InttCombinedRawDataDecoder::SetCalib (
  std::string const& domain,
  calib_load_e const& method,
  std::string const& filename)
{
  std::string str = domain;
  for(auto& c : str)
  {
    c = toupper(c);
  }
  calib_map_t::iterator itr = m_calibs.end();
  if((itr = m_calibs.find(str)) == m_calibs.end())
  {
    std::cerr << PHWHERE << "\n"
              << "\tDomain \"" << str << "\" is not settable\n"
              << "\tAvailable domains are:\n";
    ShowCalibs(std::cerr);
    std::cerr << "\t(case insensitive)" << std::endl;
    return 1;
  }

  itr->second.method = method;
  itr->second.filename = filename;

  return 0;
}

int
InttCombinedRawDataDecoder::ClearCalib (
  std::string const& domain)
{
  return SetCalib(domain, SKIP, "");
}

int
InttCombinedRawDataDecoder::ClearCalibs (
)
{
  for(auto& p : m_calibs)
  {
    p.second.method = SKIP;
    p.second.filename = "";
  }
  return 0;
}

void
InttCombinedRawDataDecoder::ShowCalibs (
  std::ostream& str)
{
  for(auto const& p : m_calibs)
  {
    str << "\t" << p.first << std::endl;
  }
}

int
InttCombinedRawDataDecoder::LoadCalibs (
)
{
  int rval = 0;
  for(auto& p : m_calibs)
  {
    int i = 0;
    switch(p.second.method)
    {
    case FROM_FILE:
      i = p.second.ptr->LoadFromFile(p.second.filename);
      break;
    case FROM_CDB:
      i = p.second.ptr->LoadFromCDB(p.second.filename);
      break;
    case SKIP:
	  if(!Verbosity())
	  {
        break;
	  }
      std::cout << PHWHERE << "\n"
                << "\tSkipping calibration \"" << p.first << "\"" << std::endl;
                << "\t(it will be used but is default-initialized)" << std::endl;
	  break;
    default:
	  break;
    }

	if(i != 0)
	{
      std::cerr << PHWHERE << "\n"
                << "\tLoading calibration \"" << p.first << "\" failed" << std::endl;
      rval = 1;
	}

	if(p.second.required && !p.second.ptr->Loaded())
	{
      std::cerr << PHWHERE << "\n"
                << "\tCalibration \"" << p.first << "\" required\n"
	            << "\tbut was not loaded" << std::endl;
      rval = 1;
	}
  }

  return rval;
}
