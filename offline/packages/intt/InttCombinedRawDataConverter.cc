#include "InttCombinedRawDataConverter.h"
#include "InttMapping.h"

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cstdlib>    // for exit
#include <iostream>   // for operator<<, endl, bas...
#include <utility>    // for pair

class PHCompositeNode;

InttCombinedRawDataConverter::InttCombinedRawDataConverter(std::string const& name)
  : SubsysReco(name)
{
  // Do nothing
}

int InttCombinedRawDataConverter::SetOutputFile(std::string const& filename)
{
  if (filename.empty())
  {
    std::cout << "InttCombinedRawDataConverter::SetOputputFile(std:: string const& filename)" << std::endl;
    std::cout << "Argument \"filename\" is empty string" << std::endl;
    std::cout << "No output was written" << std::endl;

    return 1;
  }

  std::cout << "Will write to file:" << std::endl;
  std::cout << "\t" << filename << std::endl;

  file = TFile::Open(filename.c_str(), "RECREATE");
  if (tree)
  {
    tree->SetDirectory(file);
  }

  return 0;
}

int InttCombinedRawDataConverter::WriteOutputFile()
{
  if (!file)
  {
    std::cout << "InttCombinedRawDataConverter::WriteOputputFile()" << std::endl;
    std::cout << "Member \"file\" is uninitialized" << std::endl;
    std::cout << "Did you call SetOutputFile()?" << std::endl;
    std::cout << "No output was written" << std::endl;
    return 1;
  }

  if (!tree)
  {
    std::cout << "InttCombinedRawDataConverter::WriteOputputFile()" << std::endl;
    std::cout << "Member \"tree\" is uninitialized" << std::endl;
    std::cout << "Did you call SetOutputFile()?" << std::endl;
    std::cout << "No output was written" << std::endl;
    return 1;
  }

  file->cd();
  tree->Write();
  file->Write();
  file->Close();
  return 0;
}

int InttCombinedRawDataConverter::Init(PHCompositeNode* /*topNode*/)
{
  delete tree;
  tree = new TTree("prdf_tree", "prdf_tree");
  if (file)
  {
    tree->SetDirectory(file);
  }

  tree->Branch("n_evt", &n_evt);
  tree->Branch("num_hits", &num_hits);

  branches_i =
  {
     {"flx_svr", new std::vector<Int_t>()},
     {"flx_chn", new std::vector<Int_t>()},
     {"lyr",     new std::vector<Int_t>()},
     {"ldr",     new std::vector<Int_t>()},
     {"arm",     new std::vector<Int_t>()},
     {"chp",     new std::vector<Int_t>()},
     {"chn",     new std::vector<Int_t>()},

     {"flx_bco", new std::vector<Int_t>()},
     {"adc",     new std::vector<Int_t>()},
     {"amp",     new std::vector<Int_t>()},
  };

  branches_l =
  {
     {"gtm_bco", new std::vector<Long64_t>()},
  };

  branches_d =
  {
     //{"g_x",   new std::vector<Double_t>()},
     //{"g_y",   new std::vector<Double_t>()},
     //{"g_z",   new std::vector<Double_t>()},
  };

  for (auto& itr : branches_i)tree->Branch(itr.first.c_str(), &(itr.second));
  for (auto& itr : branches_l)tree->Branch(itr.first.c_str(), &(itr.second));
  for (auto& itr : branches_d)tree->Branch(itr.first.c_str(), &(itr.second));

  if (Verbosity() > 20)
  {
    std::cout << "int InttCombinedRawDataConverter::Init(PHCompositeNode* /*topNode*/)" << std::endl;
  }
  if (Verbosity() > 20)
  {
    std::cout << "\tDone";
  }

  return 0;
}

int InttCombinedRawDataConverter::InitRun(PHCompositeNode* /*topNode*/)
{
  n_evt = 0;
  for (auto& itr : branches_i)itr.second->clear();
  for (auto& itr : branches_l)itr.second->clear();
  for (auto& itr : branches_d)itr.second->clear();

  if (Verbosity() > 20)
  {
    std::cout << "int InttCombinedRawDataConverter::InitRun(PHCompositeNode* /*topNode*/)" << std::endl;
  }
  if (Verbosity() > 20)
  {
    std::cout << "\tDone";
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataConverter::process_event(PHCompositeNode* topNode)
{
  if (!tree)
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  if (Verbosity() > 20)
  {
    std::cout << "int InttCombinedRawDataConverter::process_event(PHCompositeNode* topNode)" << std::endl;
  }
  if (Verbosity() > 20)
  {
    std::cout << "\t" << n_evt << std::endl;
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << "int InttCombinedRawDataConverter::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\t" << PHWHERE << " could not find node " << m_InttRawNodeName << std::endl;
    std::cout << "\tExiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }

  for (auto& itr : branches_i)itr.second->clear();
  for (auto& itr : branches_l)itr.second->clear();
  for (auto& itr : branches_d)itr.second->clear();

  Intt::RawData_s raw;
  Intt::Online_s onl;
  std::map<std::tuple<int, int, int, int, int>, char> hits;
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit* intthit = inttcont->get_hit(i);

    raw.felix_server = Intt::FelixFromPacket(intthit->get_packetid());
    raw.felix_channel = intthit->get_fee();
    raw.chip = intthit->get_chip_id();
    raw.channel = intthit->get_channel_id();
    onl = Intt::ToOnline(raw);

    std::tuple<int, int, int, int, int> tpl;
    std::get<0>(tpl) = onl.lyr;
    std::get<1>(tpl) = onl.ldr;
    std::get<2>(tpl) = onl.arm;
    std::get<3>(tpl) = onl.chp;
    std::get<4>(tpl) = onl.chn;

    if(hits.find(tpl) != hits.end())continue;
    hits[tpl] = 0;

    branches_i["flx_svr"]->push_back(raw.felix_server);
    branches_i["flx_chn"]->push_back(raw.felix_channel);

    branches_i["lyr"]->push_back(onl.lyr);
    branches_i["ldr"]->push_back(onl.ldr);
    branches_i["arm"]->push_back(onl.arm);
    branches_i["chp"]->push_back(onl.chp);
    branches_i["chn"]->push_back(onl.chn);

    branches_i["flx_bco"]->push_back(intthit->get_FPHX_BCO());
    branches_i["adc"]->push_back(intthit->get_adc());
    branches_i["amp"]->push_back(intthit->get_amplitude());

    branches_l["gtm_bco"]->push_back(intthit->get_bco());
  }

  num_hits = branches_l.begin()->second->size();
  tree->Fill();
  ++n_evt;

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataConverter::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
