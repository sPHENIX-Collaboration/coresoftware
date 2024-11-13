#include "InttArborist.h"

#include <Event/Event.h>
#include <Event/packet.h>

#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>

InttArborist::InttArborist()
  : m_small_branches{{"adc", new std::vector<small_t>},
                     {"amp", new std::vector<small_t>},
                     {"bco", new std::vector<small_t>},
                     {"chn", new std::vector<small_t>},
                     {"chp", new std::vector<small_t>},
                     {"fee", new std::vector<small_t>},
                     {"pid", new std::vector<small_t>}}
  , m_large_branches{{"gtm", 0}}
{
}

InttArborist::~InttArborist()
{
  for (auto& itr : m_small_branches)
  {
    delete itr.second;
  }

  delete m_tree;
}

int InttArborist::CreateOutputFile(
    std::string const& file_name)
{
  if (file_name.empty())
  {
    std::cout << "InttArborist::CreateOutputFile\n"
              << "\targument \"file_name\" is empty string" << std::endl;
    return EXIT_FAILURE;
  }

  if (m_file)
  {
    m_file->Close();
  }
  m_file = TFile::Open(file_name.c_str(), "RECREATE");
  if (!m_file)
  {
    std::cout << "InttArborist::CreateOutputFile\n"
              << "\tfailed to (re)create file " << file_name << std::endl;
    return EXIT_FAILURE;
  }

  delete m_tree;
  m_tree = new TTree(m_tree_name.c_str(), m_tree_name.c_str());
  m_file->cd();
  m_tree->SetDirectory(m_file);

  // container branches
  for (auto& itr : m_small_branches)
  {
    m_tree->Branch(itr.first.c_str(), &itr.second);
  }

  // POD-type branches
  for (auto& itr : m_large_branches)
  {
    m_tree->Branch(itr.first.c_str(), &itr.second);
  }

  return EXIT_SUCCESS;
}

int InttArborist::WriteOutputFile()
{
  if (!m_file)
  {
    std::cout << "InttArborist::WriteOutputFile\n"
              << "\tmember \"m_file\" nullptr at call" << std::endl;
    return EXIT_FAILURE;
  }

  if (!m_tree)
  {
    std::cout << "InttArborist::WriteOutputFile\n"
              << "\tmember \"m_tree\" nullptr at call" << std::endl;
    return EXIT_FAILURE;
  }

  m_file->cd();
  m_tree->Write();
  m_file->Write();
  m_file->Close();

  return EXIT_SUCCESS;
}

int InttArborist::InitRun(
    PHCompositeNode* top_node)
{
  if (!m_file)
  {
    std::cout << "InttArborist::process_event\n"
              << "\tmember \"m_file\" nullptr at call" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_tree)
  {
    std::cout << "InttArborist::process_event\n"
              << "\tmember \"m_tree\" nullptr at call" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!top_node)
  {
    std::cout << "InttArborist::process_event\n"
              << "\targument \"top_node\" nullptr at call" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttArborist::process_event(
    PHCompositeNode* top_node)
{
  if (!top_node)
  {
    std::cout << "InttArborist::process_event\n"
              << "\targument \"top_node\" nullptr at call" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  InttRawHitContainer* intt_cont = findNode::getClass<InttRawHitContainer>(
      top_node, m_intt_raw_node_name);

  if (!intt_cont)
  {
    std::cout << "InttArborist::process_event\n"
              << "\targument \"top_node\" nullptr at call" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  Gl1RawHit* gl1 = findNode::getClass<Gl1RawHit>(
      top_node, m_gl1_raw_node_name);

  m_large_branches["gtm"] = gl1 ? gl1->get_bco() : 0;
  m_large_branches["gtm"] &= (large_t{1} << 40U) - 1;

  m_clone_map.clear();
  for (auto& itr : m_small_branches)
  {
    itr.second->clear();
  }
  for (unsigned int i = 0; i < intt_cont->get_nhits(); ++i)
  {
    InttRawHit* intt_hit = intt_cont->get_hit(i);
    if (!intt_hit)
    {
      continue;
    }

    if (!m_large_branches["gtm"])
    {
      m_large_branches["gtm"] = intt_hit->get_bco();
      m_large_branches["gtm"] &= (large_t{1} << 40U) - 1;
    }

    InttMap::RawData_s raw{
        .pid = intt_hit->get_packetid(),
        .fee = intt_hit->get_fee(),
        .chp = (intt_hit->get_chip_id() + 25) % 26,
        .chn = intt_hit->get_channel_id()};

    clone_map_t::iterator clone_itr;
    if ((clone_itr = m_clone_map.find(raw)) != m_clone_map.end())
    {
      // can update our vector-valued branches at position clone_itr->second
      // for example, keep parameters of the max adc hit
      // for now I'm just doing a continue so there are no cloned hits
      // (could even make a counter for clone hits removed)
      continue;
    }

    m_small_branches["adc"]->push_back(intt_hit->get_adc());
    m_small_branches["amp"]->push_back(intt_hit->get_amplitude());
    m_small_branches["bco"]->push_back(intt_hit->get_FPHX_BCO());

    m_small_branches["chn"]->push_back(raw.chn);
    m_small_branches["chp"]->push_back(raw.chp);
    m_small_branches["fee"]->push_back(raw.fee);
    m_small_branches["pid"]->push_back(raw.pid);

    m_clone_map.insert({raw, i});
  }

  m_tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttArborist::EndRun(
    int const /*unused*/
)
{
  if (WriteOutputFile())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
