#include "MvtxFakeHitRate.h"

#include <mvtx/MvtxHitMap.h>
#include <mvtx/MvtxPixelDefs.h>
#include <mvtx/MvtxPixelMask.h>

#include <ffarawobjects/MvtxRawEvtHeader.h>
#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH1.h>
#include <TTree.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// MvtxFakeHitRate class
//==============================================================================
MvtxFakeHitRate::~MvtxFakeHitRate()
{
  // clean up
  
  
    delete m_hot_pixel_mask;
  
  
  
    delete m_hit_map;
  
}

int MvtxFakeHitRate::InitRun(PHCompositeNode* /*topNode*/)
{
  // Load the hot pixel map from the CDB
  delete m_hot_pixel_mask;  // make cppcheck happy
  m_hot_pixel_mask = new MvtxPixelMask();
  if (m_load_from_cdb)
  {
    std::cout << "MvtxFakeHitRate::InitRun - Loading hot pixel map from CDB" << std::endl;
    m_hot_pixel_mask->load_from_CDB();
  }

  // Create the hit map
  m_hit_map = new MvtxHitMap();
  m_hit_map->clear();

  std::cout << "MvtxFakeHitRate::InitRun - Writing output to " << m_outputfile << std::endl;

  // Create the output file and trees
  PHTFileServer::open(m_outputfile, "RECREATE");
  // main tree
  m_tree = new TTree("masked_pixels", "masked_pixels");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);
  m_tree->Branch("num_strobes", &m_num_strobes);
  m_tree->Branch("num_masked_pixels", &m_num_masked_pixels);
  m_tree->Branch("noise_threshold", &m_noise_threshold);
  m_tree->Branch("masked_pixels", &m_masked_pixels);

  // current mask tree
  m_current_mask = new TTree("current_mask", "current_mask");
  m_current_mask->OptimizeBaskets();
  m_current_mask->SetAutoSave(-5e6);
  m_current_mask->Branch("nmasked", &m_current_nmasked);
  m_current_mask->Branch("threshold", &m_current_threshold);
  m_current_mask->Branch("masked_pixels", &m_current_masked_pixels);

  // threshold vs nmasked histogram
  m_threshold_vs_nmasked = new TH1D("threshold_vs_nmasked", "threshold_vs_nmasked", m_max_masked_pixels, 0.5, (1.0 * m_max_masked_pixels) + 0.5);

  m_num_strobes = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxFakeHitRate::get_nodes(PHCompositeNode* topNode)
{
  // get dst nodes
  m_mvtx_raw_event_header = findNode::getClass<MvtxRawEvtHeader>(topNode, "MVTXRAWEVTHEADER");
  if (!m_mvtx_raw_event_header)
  {
    std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTXRAWEVTHEADER from Node Tree" << std::endl;
    exit(1);
  }
  if (Verbosity() > 2)
  {
    m_mvtx_raw_event_header->identify();
  }

  m_mvtx_raw_hit_container = findNode::getClass<MvtxRawHitContainer>(topNode, "MVTXRAWHIT");
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

int MvtxFakeHitRate::process_event(PHCompositeNode* topNode)
{
  // get the nodes
  if (m_num_strobes % 100000 == 0 && Verbosity() > 0)
  {
    std::cout << "MvtxFakeHitRate::process_event - processing strobe number: " << m_num_strobes << std::endl;
  }
  get_nodes(topNode);

  // get lls from the MVTX raw event header
  for (unsigned int ihit = 0; ihit < m_mvtx_raw_hit_container->get_nhits(); ihit++)
  {
    // get this hit
    auto *mvtx_hit = m_mvtx_raw_hit_container->get_hit(ihit);
    if (!mvtx_hit)
    {
      std::cout << PHWHERE << "::" << __func__ << ": Could not get MVTX hit from container. Hit index: " << ihit << std::endl;
      continue;
    }
    if (m_hot_pixel_mask->is_masked(mvtx_hit))
    {
      m_masked_pixels_in_file = true;
    }
    // get the hit info
    uint64_t strobe = mvtx_hit->get_bco();
    uint8_t layer = mvtx_hit->get_layer_id();
    if (strobe > m_last_strobe)
    {
      m_last_strobe = strobe;
      m_num_strobes++;
    }

    // if we are only looking at a specific layer, skip the rest
    if ((m_target_layer >= 0 && m_target_layer < 3) && (layer != m_target_layer))
    {
      continue;
    }

    // generate pixel key
    MvtxPixelDefs::pixelkey this_pixelkey = MvtxPixelDefs::gen_pixelkey(mvtx_hit);

    // add the hit to the hit map
    m_hit_map->add_hit(this_pixelkey, 1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxFakeHitRate::FillCurrentMaskTree()
{
  // Fill the current mask tree
  m_current_masked_pixels.clear();
  m_current_nmasked = 0;
  std::vector<MvtxPixelDefs::pixelkey> hot_pixel_map = m_hot_pixel_mask->get_hot_pixel_map();

  unsigned int nhits = m_hit_map->sum_hits(0);
  for (auto key : hot_pixel_map)
  {
    m_current_masked_pixels.push_back(key);
    m_current_nmasked++;

    unsigned int nmasked_hits = m_hit_map->get_nhits(key);
    nhits -= nmasked_hits;
  }

  m_current_threshold = calc_threshold(nhits);
  std::cout << "MvtxFakeHitRate::FillCurrentMaskTree - nmasked: " << m_current_nmasked << " threshold: " << m_current_threshold << std::endl;

  m_current_mask->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

double MvtxFakeHitRate::calc_threshold(int nhits) const
{
  // Calculate the noise threshold
  if (nhits == 0)
  {
    return 0.0;
  }
  double npixels = 226492416.0;
  double denom = npixels * m_num_strobes;
  if (denom == 0)
  {
    return 0.0;
  }
  return (static_cast<double>(nhits) / denom);
}

int MvtxFakeHitRate::CalcFHR()
{
  // Calculate the fake hit rate

  MvtxHitMap::pixel_hit_vector_t pixel_hit_vector = m_hit_map->get_pixel_hit_vector();

  int npixels_with_hits = pixel_hit_vector.size();
  if (npixels_with_hits == 0)
  {
    std::cout << "MvtxFakeHitRate::CalcFHR - No pixels with hits" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_max_masked_pixels = std::min(npixels_with_hits, m_max_masked_pixels);

  // sort the pixel hit vector by hit count
  std::sort(pixel_hit_vector.begin(), pixel_hit_vector.end(), [](const MvtxHitMap::pixel_hits_pair_t& a, const MvtxHitMap::pixel_hits_pair_t& b)
            { return a.second > b.second; });

  // get initial values
  unsigned int nhits_all = 0;
  m_masked_pixels.clear();
  for (auto & it : pixel_hit_vector)
  {
    nhits_all += it.second;
  }
  m_noise_threshold = calc_threshold(nhits_all);
  m_threshold_vs_nmasked->Fill(0.5, m_noise_threshold);
  m_tree->Fill();

  for (int nmasked = 0; nmasked < m_max_masked_pixels; nmasked++)
  {
    m_masked_pixels.clear();
    m_num_masked_pixels = 0;
    int ipixel = 0;
    int nhits = nhits_all;

    for (auto & it : pixel_hit_vector)
    {
      if (ipixel < nmasked)
      {
        m_masked_pixels.push_back(it.first);
        m_num_masked_pixels++;
        nhits -= it.second;
      }
      ipixel++;
    }
    if (nhits < 0)
    {
      std::cout << "MvtxFakeHitRate::CalcFHR - Error: nhits < 0" << std::endl;
      nhits = 0;
    }
    m_noise_threshold = calc_threshold(nhits);

    m_threshold_vs_nmasked->Fill(m_num_masked_pixels + 0.5, m_noise_threshold);
    if (m_num_masked_pixels % 1000 == 0)
    {
      std::cout << "MvtxFakeHitRate::CalcFHR - nmasked: " << m_num_masked_pixels << " threshold: " << m_noise_threshold << std::endl;
    }
    m_tree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxFakeHitRate::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "MvtxFakeHitRate::End - Number of strobes: " << m_num_strobes << std::endl;
  // Fill the current mask tree
  FillCurrentMaskTree();
  if (m_masked_pixels_in_file)
  {
    std::cout << "MvtxFakeHitRate::End - Masked pixels found in file" << std::endl;
  }
  else
  {
    std::cout << "MvtxFakeHitRate::End - No masked pixels found in file" << std::endl;
  }

  // Calculate the fake hit rate
  CalcFHR();

  // Write the output
  std::cout << "MvtxFakeHitRate::End - Writing output to " << m_outputfile << std::endl;
  PHTFileServer::cd(m_outputfile);
  m_tree->Write();
  m_current_mask->Write();
  m_threshold_vs_nmasked->GetXaxis()->SetNdivisions(505);
  m_threshold_vs_nmasked->GetYaxis()->SetNdivisions(505);
  m_threshold_vs_nmasked->GetXaxis()->SetTitle("Masked pixels");
  m_threshold_vs_nmasked->GetYaxis()->SetTitle("(Active pix. / tot. pix.) per strobe");
  m_threshold_vs_nmasked->SetFillColor(kRed);
  m_threshold_vs_nmasked->SetLineColor(kRed);
  m_threshold_vs_nmasked->SetMarkerColor(kRed);
  m_threshold_vs_nmasked->Draw("PE1");
  m_threshold_vs_nmasked->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}
