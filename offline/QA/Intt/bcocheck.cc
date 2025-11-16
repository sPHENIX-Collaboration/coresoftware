#include "bcocheck.h"

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cstdint>
#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>

bcocheck::bcocheck(const std::string &name, const int run_num, const int felix_num)
  : SubsysReco(name)
  , run_num_(run_num)
  , felix_num_(felix_num)
{
  std::cout << "felix_num_=" << felix_num_ << "felix_num=" << felix_num << std::endl;
}

int bcocheck::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in bcocheck" << std::endl;
  }

  return 0;
}

int bcocheck::InitRun(PHCompositeNode *topNode)
{
  if (!topNode)
  {
    std::cout << "bcocheck::InitRun(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;

    return 1;
  }

  std::cout << "felix_num_=" << felix_num_ << std::endl;

  // Initialize histograms
  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    if (felix != felix_num_)
    {
      continue;
    }

    std::string const name = "h2_bco_felix_" + std::to_string(felix);
    std::string const title = name + std::format("_Run{}", run_num_);
    h_full_bco[felix] = new TH1D(name.c_str(), title.c_str(), 128, 0, 128);
    h_full_bco[felix]->SetXTitle("BCO_FULL - BCO");
    h_full_bco[felix]->SetMinimum(0);
    tf_output_[felix] = new TFile(std::format("./bco_000{}_intt{}.root", run_num_, felix).c_str(), "RECREATE");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int bcocheck::process_event(PHCompositeNode *topNode)
{
  // std::cout<<"1234"<<std::endl;

  // get raw hit
  std::string const m_InttRawNodeName = "INTTRAWHIT";
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "bcocheck::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  uint64_t const longbco_full = inttcont->get_nhits() > 0
                                    ? inttcont->get_hit(0)->get_bco()
                                    : std::numeric_limits<uint64_t>::max();
  // uint64_t difevent_bcofull = (longbco_full &bit )-(long_prev_bcofull &bit);
  // h_interval->Fill(difevent_bcofull);
  uint64_t const bco_full = longbco_full & 0x7FU;
  // std::cout<<"longbco_full="<<longbco_full<<std::endl;

  ievent_++;
  // std::cout<<"5678"<<std::endl;
  /*if(ievent_<500000){
    return Fun4AllReturnCodes::EVENT_OK;
  }*/

  if ((ievent_ % 100 == 0 && ievent_ < 1000) || ievent_ % 1000 == 0)
  {
    std::cout << "Process event #" << ievent_ << std::endl;
  }

  int const nhits = inttcont->get_nhits();

  for (int i = 0; i < nhits; i++)
  {
    InttRawHit *intthit = inttcont->get_hit(i);

    int const fnum = intthit->get_packetid() - 3001;  // packet id
    int const bco = intthit->get_FPHX_BCO();          // FPHX bco

    int bco_diff = bco_full - bco;
    if (bco_diff < 0)
    {
      bco_diff = bco_diff + 128;
    }

    h_full_bco[fnum]->Fill(bco_diff);
  }
  // std::cout<<"910"<<std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int bcocheck::End(PHCompositeNode * /*topNode*/)
{
  // Mixupfraction();

  if (Verbosity() > 1)
  {
    std::cout << "Processing bcocheck done" << std::endl;
  }

  this->DrawHists();

  return Fun4AllReturnCodes::EVENT_OK;
}

void bcocheck::DrawHists()
{
  gStyle->SetOptStat(0);

  for (int felix = 0; felix < kFelix_num_; felix++)
  {
    if (felix != felix_num_)
    {
      continue;
    }

    tf_output_[felix]->cd();

    h_full_bco[felix]->Write();

    tf_output_[felix]->Close();
  }
}
