#include "MvtxHitMapWriter.h"

#include <mvtx/MvtxHitMap.h>
#include <mvtx/MvtxPixelDefs.h>
#include <mvtx/MvtxPixelMask.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <ffarawobjects/MvtxRawEvtHeader.h>
#include <ffarawobjects/MvtxRawEvtHeaderv2.h>
#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>
#include <ffarawobjects/MvtxRawHitContainerv1.h>
#include <ffarawobjects/MvtxRawHitv1.h>

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <TTree.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// MvtxHitMapWriter class
//==============================================================================
int MvtxHitMapWriter::InitRun(PHCompositeNode* /*topNode*/)
{
  // Create the hit map
  delete m_hit_map;  // make cppcheck happy
  m_hit_map = new MvtxHitMap();
  m_hit_map->clear();

  std::cout << "MvtxHitMapWriter::InitRun - Writing output to " << m_outputfile << std::endl;

  // Create the output file and trees
  PHTFileServer::get().open(m_outputfile, "RECREATE");
  // main tree
  m_tree_info = new TTree("info", "MVTX Hit Map Info");
  m_tree_info->Branch("num_strobes", &m_num_strobes, "num_strobes/i");
  m_tree_info->Branch("nhits_total", &m_nhits_total, "nhits_total/i");
  m_tree_info->Branch("num_fired_pixels", &m_num_fired_pixels, "num_fired_pixels/i");

  // hit map tree
  m_tree = new TTree("hit_map", "MVTX Hit Map");
  m_tree->Branch("pixels", &m_pixels);
  m_tree->Branch("nhits", &m_nhits);

  m_num_strobes = 0;
  m_nhits_total = 0;
  m_num_fired_pixels = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxHitMapWriter::get_nodes(PHCompositeNode* topNode)
{
  // get dst nodes
  m_mvtx_raw_event_header = findNode::getClass<MvtxRawEvtHeaderv2>(topNode, "MVTXRAWEVTHEADER");
  if (!m_mvtx_raw_event_header)
  {
    std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTXRAWEVTHEADER from Node Tree" << std::endl;
    exit(1);
  }
  if (Verbosity() > 2)
  {
    m_mvtx_raw_event_header->identify();
  }

  m_mvtx_raw_hit_container = findNode::getClass<MvtxRawHitContainerv1>(topNode, "MVTXRAWHIT");
  if (!m_mvtx_raw_hit_container)
  {
    std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTXRAWHIT from Node Tree" << std::endl;
    exit(1);
  }
  if (Verbosity() > 2)
  {
    m_mvtx_raw_hit_container->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxHitMapWriter::process_event(PHCompositeNode* topNode)
{
  // get the nodes
  if (m_num_strobes % 100000 == 0)
  {
    std::cout << "MvtxHitMapWriter::process_event - processing strobe number: " << m_num_strobes << std::endl;
  }

  get_nodes(topNode);

  // get lls from the MVTX raw event header
  for (unsigned int ihit = 0; ihit < m_mvtx_raw_hit_container->get_nhits(); ihit++)
  {
    // get this hit
    auto mvtx_hit = m_mvtx_raw_hit_container->get_hit(ihit);
    if (!mvtx_hit)
    {
      std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTX hit from container. Hit index: " << ihit << std::endl;
      continue;
    }

    // get the hit info
    uint64_t strobe = mvtx_hit->get_bco();
    if (strobe > m_last_strobe)
    {
      m_last_strobe = strobe;
      m_num_strobes++;
    }

    // generate pixel key
    MvtxPixelDefs::pixelkey this_pixelkey = MvtxPixelDefs::gen_pixelkey(mvtx_hit);

    // add the hit to the hit map
    m_hit_map->add_hit(this_pixelkey, 1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxHitMapWriter::FillTree()
{
  // get the hit map
  MvtxHitMap::pixel_hit_vector_t pixel_hit_vector = m_hit_map->get_pixel_hit_vector();
  std::sort(pixel_hit_vector.begin(), pixel_hit_vector.end(), [](const MvtxHitMap::pixel_hits_pair_t& a, const MvtxHitMap::pixel_hits_pair_t& b)
            { return a.second > b.second; });

  m_num_fired_pixels = pixel_hit_vector.size();
  if (m_num_fired_pixels == 0)
  {
    std::cout << "MvtxHitMapWriter::FillTree - No pixels with hits" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_nhits_total = 0;
  m_pixels.clear();
  m_nhits.clear();

  for (auto it = pixel_hit_vector.begin(); it != pixel_hit_vector.end(); ++it)
  {
    m_pixels.push_back(it->first);
    m_nhits.push_back(it->second);
    m_nhits_total += it->second;
  }

  m_tree->Fill();
  m_tree_info->Fill();

  std::cout << "MvtxHitMapWriter::FillTree - Number of pixels with hits: " << m_num_fired_pixels << std::endl;
  std::cout << "MvtxHitMapWriter::FillTree - Total number of hits: " << m_nhits_total << std::endl;
  std::cout << "MvtxHitMapWriter::FillTree - Number of strobes: " << m_num_strobes << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxHitMapWriter::End(PHCompositeNode* /*topNode*/)
{
  // Fill the current mask tree
  FillTree();
  std::cout << "MvtxHitMapWriter::End - Writing output to " << m_outputfile << std::endl;
  PHTFileServer::get().cd(m_outputfile);
  PHTFileServer::get().write(m_outputfile);
  return Fun4AllReturnCodes::EVENT_OK;
}
